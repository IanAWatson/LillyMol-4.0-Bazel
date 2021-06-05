# frozen_string_literal: true

# iwdescr fingerprint class for gfp_make

require_relative "lib/fp_common.rb"

class DSC
  @@rx = Regexp.new("^DSC")
  @@description = "descriptor based fingerprint"
  @@executable = "iwdescr"

  def description
    return @@description
  end

  def match?(fp)
    return @@rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:)
    m = /^DSC(\d+)*/.match(fp)
    if not m
      raise "Unrecognized DSC fp form '#{fp}'"
    end

    cmd = FpCommon.initial_command_stem(@@executable, first_in_pipeline:first_in_pipeline,
                extra_qualifiers:extra_qualifiers)
    cmd
  end
  
end
