# frozen_string_literal: true

# MaccsKeys class for gfp_make

require_relative 'lib/fp_common.rb'

class MK
  @@rx = Regexp.new("^MK")
  @@description = "MACCS keys fingerprint"
  @@executable = "maccskeys_fn5"

  def description
    return @@description
  end

  def match?(fp)
    return @@rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:)
    m = /^MK(\d+)*/.match(fp)
    if not m
      raise "Unrecognized MK fp form '#{fp}'"
    end

    cmd = FpCommon.initial_command_stem(@@executable, first_in_pipeline:first_in_pipeline,
                extra_qualifiers:extra_qualifiers)
    cmd
  end
end
