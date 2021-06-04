# frozen_string_literal: true

# IWFP class for gfp_make

require_relative "lib/fp_common.rb"

class IW
  @@rx = Regexp.new("^IW")
  @@description = "Linear path fingerprint"
  @@executable = "iwfp"

  def description
    return @@description
  end

  def match?(fp)
    return @@rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:)
    m = /^IW(\d+)*/.match(fp)
    if not m
      raise "Unrecognized IW fp form '#{fp}'"
    end

    cmd = FpCommon.initial_command_stem(@@executable, first_in_pipeline:first_in_pipeline,
                extra_qualifiers:extra_qualifiers)
    path_length, atype = FpCommon.parse_fp_token(fp[2..])

    cmd << ' -J FPIW'
    cmd << "#{path_length} -R #{path_length}" if path_length
    cmd << " -P #{atype}" if atype
    cmd
  end
  
end
