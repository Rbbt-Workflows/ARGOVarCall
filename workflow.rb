require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/ARGOVarCall'

Workflow.require_workflow "HTS"
module ARGOVarCall
  module CompareIndels
    extend Workflow

    input :combined_VCF, :file, "Combined caller VCF", nil, :nofile => true
    task :mutation_positions => :tsv do |vcf|

      parser = TSV::Parser.new vcf, :type => :list
      dumper = TSV::Dumper.new :key_field => "Mutation ID", :fields => ["Start", "Delete", "Added", "Genomic mutation", "CHR"] + parser.fields, :type => :list
      dumper.init
      TSV.traverse parser, :type => :array, :into => dumper do |chr,values|
        parts = [chr] + values
        chr, pos, rsid, ref, alts, qual, filter = parts
        pos = pos.to_i

        filter = filter.split(";").collect{|v| v.split("--").first } * "+"

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
      TSV.traverse step(:mutation_positions), :into => :stream do |id,values|
        chr = id.split(":").first
        pos, ins, del = values[0..2].collect{|v| v.to_i }
        eend = pos + ins + del
        [chr, pos, ins, del, id] * "\t"
      end
    end

    dep :ranges
    task :overlaps => :tsv do
      io = CMD.cmd('sort -k1,1 -k2,2n', :in => step(:ranges).join.path.open, :pipe => true)

      chunk = []
      chunk_start, chunk_end, chunk_chr = nil, nil, nil
      slack = 1
      tsv = TSV.setup({}, "Range~Mutation ID#:type=:flat")
      TSV.traverse io, :type => :array do |line|
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

    def self.dna_after_mutations(chr_sequence, mutations)
      chr_sequence = chr_sequence.dup.split("")
      changes = mutations.collect{|m| parts = m.split(":"); [parts[1].to_i - 1, parts[2]]}
      changes.sort_by{|p| p.first }.reverse.each do |pos,change|
        change = change.split("")

        if change[0] == "-"
          while change[0] == "-"
            chr_sequence.delete_at pos
            change.shift
          end
          if Array === chr_sequence[pos]
            Log.warn "Double hit: #{mutations.inspect}"
            next
          end
          chr_sequence[pos] = change * "" + chr_sequence[pos] if change.any?
        else
          change[0] = chr_sequence[pos] if change[0] == "+"
          chr_sequence[pos] = change if change.any?
        end

      end

      chr_sequence = chr_sequence*""

      chr_sequence.gsub("-",'')
    end

    helper :chromosome do |code|
      @chromosomes ||= {}
      @chromosomes[code] ||= begin
                               require 'rbbt/sources/organism'
                               organism = "Hsa/may2017"
                               Organism.send("chromosome_#{code.sub('chr','')}", organism).read
                             end
    end


    dep :overlaps
    task :reconstructed_sequence => :tsv do
      pad = 3
      range_sequences = {}
      overlaps = step(:overlaps).path.tsv

      TSV.traverse overlaps, :bar => self.progress_bar("Processing overlaps") do |range, values|
        range = range.first if Array === range

        chr, start, eend = range.split(":")
        start = start.to_i
        eend = eend.to_i

        chr_sequence = begin
                         chromosome(chr).dup
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
          vc_sequences[vc] = ARGOVarCall.dna_after_mutations(seq, changes)
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
    task :extended_validations => :tsv do
      result = TSV.setup({}, :key_fields => "Range", :fields => %w(Mutation IDs), :type => :flat) 
      TSV.traverse step(:reconstructed_sequence) do |range,values,fields|
        mutations, reference, caller_sequences = values

        seq_vcs = {}
        fields.zip(values)[3..-1].each do |vc,seq|
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

  end

  extend Workflow

  input :vcf_list, :array, "Array of VCF files to process", nil, :required => true
  input :caller_names, :array, "Corresponding names for the files (must have same length)"
  extension :vcf
  task :combine_vcf => :text do |vcf_list, caller_names|
    if caller_names
      list = Misc.zip2hash(caller_names, vcf_list)
    else
      list = Misc.zip2hash(vcf_list.collect{|f| File.basename(f).sub(/\.vcf(\.gz)?$/i,'')}, vcf_list)
    end

    HTS.combine_caller_vcfs(list)
  end


  dep_task :validations, CompareIndels, :extended_validations
end

  #require 'ARGOVarCall/tasks/basic.rb'

  #require 'rbbt/knowledge_base/ARGOVarCall'
  #require 'rbbt/entity/ARGOVarCall'

