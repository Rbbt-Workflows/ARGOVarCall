module ARGOVarCall
  def self.filter_overlaps(mutations)
    filtered = []
    eend = nil
    mutation_info = mutations.collect do |m| 
      p = m.split(":")
      size = p.last.scan(/-*/).first.length
      p + [size]
    end

    last = nil
    mutation_info.sort{|p| - p[1].to_i }.reverse.collect do |p|
      chr, pos, alt, size = p
      next if last && pos.to_i + size.to_i - 1 >= last
      last = pos.to_i
      [chr, pos, alt] * ":"
    end.compact
  end

  def self.dna_after_mutations(chr_sequence, mutations)
    chr_sequence = chr_sequence.dup.split("")
    changes = filter_overlaps(mutations).collect{|m| parts = m.split(":"); [parts[1].to_i - 1, parts[2]]}
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
      chr_sequence = (chr_sequence*"").split("")

    end

    chr_sequence = chr_sequence*""

    chr_sequence.gsub("-",'')
  end

end
