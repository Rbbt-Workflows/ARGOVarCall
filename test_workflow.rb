require 'rbbt-util'
require 'rbbt/workflow'
require 'test/unit'

class TestWorkflow < Test::Unit::TestCase
  def workflow
    @@workflow ||= Workflow.require_workflow __FILE__.sub('test_','')
  end

  def last_job
    task_name = workflow.tasks.keys.last
    workflow.job(task_name)
  end

  def first_job
    task_name = workflow.tasks.keys.first
    workflow.job(task_name)
  end

  def test_true
    workflow
    chr_seq = "123456789"
    mutations = "chr1:7:-"

    assert_equal "12345689", ARGOVarCall.dna_after_mutations(chr_seq, ["chr1:7:-"])
    assert_equal "123456", ARGOVarCall.dna_after_mutations(chr_seq, ["chr1:7:---"])
    assert_equal "123456A89", ARGOVarCall.dna_after_mutations(chr_seq, ["chr1:7:A"])
    assert_equal "1234567A89", ARGOVarCall.dna_after_mutations(chr_seq, ["chr1:7:+A"])
    assert_equal "123456ABC89", ARGOVarCall.dna_after_mutations(chr_seq, ["chr1:7:ABC"])
    assert_equal "123456AB89", ARGOVarCall.dna_after_mutations(chr_seq, ["chr1:7:-AB"])
    assert_equal "123456AB9", ARGOVarCall.dna_after_mutations(chr_seq, ["chr1:7:--AB"])
  end
end

