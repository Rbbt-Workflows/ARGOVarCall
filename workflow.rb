require 'rbbt-util'
Misc.add_libdir if __FILE__ == $0

require 'VariantConsensus'

module VariantConsensus
  module CompareIndels
    extend Workflow

    input :combined_vcf, :file, "Combined caller VCF", nil, :nofile => true
    task :mutation_positions => :tsv do |vcf|

      parser = TSV::Parser.new vcf, :type => :list
      dumper = TSV::Dumper.new :key_field => "Mutation ID", :fields => ["Start", "Delete", "Added", "Genomic mutation", "CHR"] + parser.fields, :type => :list
      dumper.init
      TSV.traverse parser, :into => dumper, :bar => self.progress_bar("Finding mutation positions") do |chr,values|
        parts = [chr] + values
        chr, pos, rsid, ref, alts, qual, filter = parts
        pos = pos.to_i

        rsid = '.'

        filter = filter.split(";").collect{|v| v.split("--").first }.uniq * "+"

        new_pos, muts = Misc.correct_vcf_mutation(pos, ref, alts)
        res = muts.collect do |mut|
          id = [chr, pos, rsid, ref, alts,new_pos,mut, filter] * ":"
          del = mut.split("").select{|b| b == "-"}.length
          ins = mut.split("").reject{|b| b == "-" || b == "+"}.length
          [id, [new_pos, ins, del, [chr,new_pos,mut] * ":", chr] + parts]
        end.compact

        res.extend MultipleResult
        res
      end
    end

    dep :mutation_positions
    task :ranges => :array do
      TSV.traverse step(:mutation_positions), :into => :stream, :bar => self.progress_bar("Calulating ranges") do |id,values|
        chr = id.split(":").first
        pos, ins, del = values[0..2].collect{|v| v.to_i }
        eend = pos + ins + del
        [chr, pos, ins, del, id] * "\t"
      end
    end

    dep :ranges
    task :overlaps => :tsv do
      step(:ranges).join
      s = step(:ranges).path.open
      io = CMD.cmd('sort -k1,1 -k2,2n', :in => s, :pipe => true)

      chunk = []
      chunk_start, chunk_end, chunk_chr = nil, nil, nil
      slack = 0
      tsv = TSV.setup({}, "Range~Mutation ID#:type=:flat")
      TSV.traverse io, :type => :array, :bar => self.progress_bar("Finding overlaps") do |line|
        chr, start, ins, del, id = line.split("\t")
        start = start.to_i
        ins = ins.to_i
        del = del.to_i

        info = [start, ins, del, id]

        if chunk.empty?
          chunk = [info]
          chunk_start = start
          chunk_end = start + del + 1
          chunk_chr = chr
          next
        end

        if chr == chunk_chr && start < chunk_end + slack
          chunk << info
          new_chunk_end = start + del + 1
          chunk_end = new_chunk_end if new_chunk_end > chunk_end
        else
          if chunk.length > 1
            chunk.first.last.split(":").first
            range = [chunk_chr,chunk_start,chunk_end] * ":"
            tsv[range] = chunk.collect{|i| i.last }
          end
          chunk = [info]
          chunk_start = start
          chunk_end = start + del + 1
          chunk_chr = chr
        end
      end

      tsv
    end

    helper :chromosome do |code,organism="Hsa/may2017"|
      @chromosomes ||= {}
      @chromosomes[code] ||= begin
                               require 'rbbt/sources/organism'
                               code = "MT" if code == "chrM"
                               Organism.send("chromosome_#{code.sub('chr','')}", organism).read
                             end
    end


    input :organism, :select, "Organism code", "Hsa/may2017"
    dep :overlaps, :compute => :produce
    task :reconstructed_sequence => :tsv do |organism|
      pad = 3
      range_sequences = {}
      overlaps = step(:overlaps).path.tsv

      TSV.traverse overlaps, :bar => self.progress_bar("Processing overlaps") do |range, values|
        range = range.first if Array === range

        chr, start, eend = range.split(":")
        next if chr.include? "_"
        start = start.to_i
        eend = eend.to_i

        chr_sequence = begin
                         chromosome(chr, organism).dup
                       rescue 
                         next
                       end

        seq_start = start - pad - 1 
        seq_end = eend + pad - 1
        seq = chr_sequence[seq_start..seq_end]

        vc_mutations = {}
        values.collect do |value|
          value.split(":").last.split("+").each do |vc|
            vc_mutations[vc] ||= []
            organism = "Hsa/feb2014"
            vc_mutations[vc] << value
          end
        end

        vc_sequences = {"Reference" => seq}
        vc_mutations.collect do |vc,mutations|
          changes = mutations.collect do |info|
            parts = info.split(":")
            chr, pos, alt = parts.values_at 0,5,6
            [chr, pos.to_i - seq_start, alt] * ":"
          end
          vc_sequences[vc] = VariantConsensus.dna_after_mutations(seq, changes)
        end
        range_sequences[range] = vc_sequences
      end

      variant_callers = range_sequences.values.collect{|h| h.keys }.flatten.uniq.sort

      result = TSV.setup({}, :key_field => "Range", :fields => ["Mutations"] + variant_callers, :type => :list)

      range_sequences.each do |range,info|
        mutations = overlaps[range]

        result[range] = [mutations] + variant_callers.collect{|vc| range_sequences[range][vc] }
      end

      result
    end

    dep :reconstructed_sequence
    task :equivalent_changes => :tsv do
      result = TSV.setup({}, :key_fields => "Range", :fields => ["Mutation IDs"], :type => :flat) 
      TSV.traverse step(:reconstructed_sequence) do |range,values,fields|
        mutations, reference, caller_sequences = values

        seq_vcs = {}
        fields.zip(values)[2..-1].each do |vc,seq|
          next if seq.empty?
          seq_vcs[seq] ||= []
          seq_vcs[seq] << vc
        end

        seq_vcs.select{|s,vcs| vcs.length > 1}.each do |s,vcs|
          good_mutations = mutations.split("|").select{|muts| (vcs & muts.split(":").last.split("+")).any? }
          good_mutations = good_mutations.sort_by{|m| m.split(":")[-2].length }
          next unless good_mutations.length > 1
          result[range] ||= []
          result[range] << vcs*"+" + " => " + good_mutations *";"
        end
      end
      result
    end

    dep :equivalent_changes
    task :equivalent_mutations => :tsv do
      dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => %w(Callers ID Others), :type => :double
      dumper.init
      TSV.traverse step(:equivalent_changes), :into => :dumper do |range,entries|
        res = []
        res.extend MultipleResult
        entries.each do |entry|
          callers, mutations = entry.split(" => ")
          id = Misc.digest(mutations)
          mutations = mutations.split(";").collect{|m| m.split(":")[0..4] * ":" }.uniq
          mutations.each do |mutation|
            res << [mutation, [callers.split("+"), id, mutations - [mutation]]]
          end
        end
        res
      end
    end

  end

  extend Workflow

  input :vcf_list, :file_array, "Array of VCF files to process", nil, :required => true
  input :caller_names, :array, "Corresponding names for the files (must have same length)"
  extension :vcf
  task :combine_vcf => :text do |vcf_list, caller_names|
    if caller_names
      list = Misc.zip2hash(caller_names, vcf_list)
    else
      list = Misc.zip2hash(vcf_list.collect{|f| File.basename(f).sub(/\.vcf(\.gz)?$/i,'')}, vcf_list)
    end

    VariantConsensus.combine_caller_vcfs(list)
  end

  dep :combine_vcf
  dep_task :equivalent_mutations, CompareIndels, :equivalent_mutations, :combined_vcf => :combine_vcf

  dep :combine_vcf
  dep :equivalent_mutations
  extension :vcf
  task :extended_validation_combined_vcf => :text do
    equivalent_mutations = step(:equivalent_mutations).load

    TSV.traverse step(:combine_vcf), :type => :line, :bar => self.progress_bar("Adding validation information"), :into => :stream do |line|
      next line if line =~ /^##/
      if line =~ /^#?CHR/
        '##INFO=<ID=UniqueRegionID,Number=1,Type=String,Description="Unique ID of region harboring consecutive or overlapping mutations, used to avoid double counting">' + "\n" +
        '##INFO=<ID=ValidatedBy,Number=.,Type=String,Description="Callers that validate the mutation, posibly with different mutations">' + "\n" +
        '##INFO=<ID=NumCallers,Number=1,Type=Integer,Description="Number of callers that call the mutation">' + "\n" +
        '##INFO=<ID=NumValidators,Number=1,Type=Integer,Description="Number of callers that validate the mutation, posibly with different mutations">' + "\n" +
          line
      else
        parts = line.split("\t")
        mutation = parts[0..4] * ":"
        info = Hash[parts[7].split(";").collect{|p| p.split("=") }]
        if equivalent_mutations.include? mutation
          callers, id, others = equivalent_mutations[mutation]
          info["ValidatedBy"] = (info["CalledBy"].split(",") + callers).uniq * ","
          info["UniqueRegionID"] = "Ambiguous-" + id.first
        else
          info["ValidatedBy"] = info["CalledBy"]
          info["UniqueRegionID"] = Misc.digest(mutation)
        end
        info["NumCallers"] = info["CalledBy"].split(",").length
        info["NumValidators"] = info["ValidatedBy"].split(",").length
        parts[7] = info.collect{|p| p * "=" } * ";"

        parts * "\t"
      end
    end
  end

  dep :extended_validation_combined_vcf
  input :consensus_field, :select, "Consensus field to use", "ValidatedBy", :select_options => %w(CalledBy ValidatedBy)
  input :min_validation, :integer, "Minimum number of validating callers", 2
  input :filter_regime, :select, "How many callers must have a PASS filter originally", :one, :select_options => %w(all one none)
  extension :vcf
  task :fully_annotated_consensus_vcf => :text do |consensus_field, min_validation,filter_regime|
    TSV.traverse step(:extended_validation_combined_vcf), :type => :line, :bar => self.progress_bar("Calculating consensus"), :into => :stream do |line|
      next line if line =~ /^##/

      if line =~ /^#?CHR/
        "##FILTER=<ID=PASS,Description=\"Consensus validated, at least #{min_validation} callers support the variant (#{filter_regime} need to have an original PASS filter)}\">" + "\n" +
        "##FILTER=<ID=FAIL,Description=\"Consensus not validated, less than #{min_validation} callers support the variant (#{filter_regime} need to have an original PASS filter)}\">" + "\n" +
        line
      else
        parts = line.split("\t")

        info = Hash[parts[7].split(";").collect{|p| p.split("=") }]

        callers = info[consensus_field].split(",")

        filters = parts[6].split(";").collect{|p| p.split("--").last }

        case filter_regime.to_s
        when 'all'
          pass = filters.uniq == ["PASS"]
        when 'one'
          pass = filters.include? "PASS"
        else
          pass = true
        end


        if pass && callers.length >= min_validation
          parts[6] = "PASS;" + parts[6]
        else
          parts[6] = "FAIL;" + parts[6]
        end

        parts * "\t"
      end
    end
  end


  dep :fully_annotated_consensus_vcf
  extension :vcf
  task :consensus_vcf => :text do
    TSV.traverse step(:fully_annotated_consensus_vcf), :type => :line, :bar => self.progress_bar("Cleaning VCF of extended annotations"), :into => :stream do |line|
      if line.start_with? "##"
        if line.include? "--"
          next
        else
          next line
        end
      end

      next line if line.start_with? "#"

      parts = line.split("\t")

      parts[6] = parts[6].split(";").reject{|p| p.include? "--" } * ";"

      info = Hash[parts[7].split(";").collect{|p| p.split("=") }]
      info.delete_if{|k,v| k.include? "--" }
      parts[7] = info.collect{|p| p * "=" } * ";"

      good_format = parts[8].split(":").collect{|p| ! p.include? "--" }
      parts[8] = good_format.zip(parts[8].split(":")).select{|g,p| g }.collect{|g,p| p } * ":"
      parts[9] = good_format.zip(parts[9].split(":")).select{|g,p| g }.collect{|g,p| p } * ":"
      parts[10] = good_format.zip(parts[10].split(":")).select{|g,p| g }.collect{|g,p| p } * ":"


      parts * "\t"
    end

  end
  
end

#require 'VariantConsensus/tasks/basic.rb'

#require 'rbbt/knowledge_base/VariantConsensus'
#require 'rbbt/entity/VariantConsensus'

