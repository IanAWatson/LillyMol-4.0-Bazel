# frozen_string_literal: true

# EC fingerprint module for gfp_make

require_relative "lib/fp_common.rb"

class EC
  @@rx = Regexp.new("^EC")
  @@description = "Circular fingerprints"
  @@executable = "iwecfp"

  def match?(fp)
    return @@rx.match?(fp)
  end

  def description
    return @@description
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:)
    m = /^EC([A-Z]*)(\d+)*:*(\S+)/.match(fp)
    if not m
      raise "Unrecognized EC fp form '#{fp}'"
    end
    cmd = FpCommon.initial_command_stem(@@executable, first_in_pipeline:first_in_pipeline,
                extra_qualifiers:extra_qualifiers)
    radius, atype = FpCommon.parse_fp_token(fp[2..])
#   $stderr << "Radius #{radius} predefined '#{predefined}' ust_atype #{ust_atype}\n"

    radius = '3' unless radius

    cmd << " -J NCEC#{radius} -R #{radius}"

    cmd << " -P #{atype}" if atype

    cmd
  end
  
end
