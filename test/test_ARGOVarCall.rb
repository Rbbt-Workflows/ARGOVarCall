require File.expand_path(__FILE__).sub(%r(/test/.*), '/test/test_helper.rb')
require File.expand_path(__FILE__).sub(%r(.*/test/), '').sub(/test_(.*)\.rb/,'\1')

require 'rbbt-util'

class TestARGOVarCall < Test::Unit::TestCase
  def test_dna_after_mutations
    seq = 'tcgatcgatcga'

    assert_equal "a" + seq[1..-1], ARGOVarCall.dna_after_mutations(seq, ['chr:1:a'])
    assert_equal "aa" + seq[2..-1], ARGOVarCall.dna_after_mutations(seq, ['chr:1:--aa'])
    assert_equal "aa" + seq[2..-1], ARGOVarCall.dna_after_mutations(seq, ['chr:1:a', 'chr:2:a'])

    assert_equal "aa" + seq[2..-1], ARGOVarCall.dna_after_mutations(seq, ['chr:1:a', 'chr:2:a', 'chr:1:--aa'])
    assert_equal "aa" + seq[2..-1], ARGOVarCall.dna_after_mutations(seq, ['chr:1:a', 'chr:2:a', 'chr:1:--aa'])

    assert_equal "c" + seq[3..-1], ARGOVarCall.dna_after_mutations(seq, ['chr:1:--', 'chr:3:c'])
    assert_equal "c" + seq[3..-1], ARGOVarCall.dna_after_mutations(seq, ['chr:1:---c'])
    assert_equal "c" + seq[3..-1], ARGOVarCall.dna_after_mutations(seq, ['chr:1:--', 'chr:3:c', 'chr:1:---c'])
  end
end

