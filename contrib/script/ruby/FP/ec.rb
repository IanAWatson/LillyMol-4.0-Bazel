# frozen_string_literal: true

# EC fingerprint module for gfp_make

require_relative 'lib/fp_common'

# EC fingerprint.
class EC
  attr_reader :description

  def initialize
    @rx = Regexp.new('^EC')
    @description = 'Circular fingerprints'
    @executable = 'iwecfp'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^EC([A-Z]*)(\d+)*:*(\S+)*/.match(fp)
    raise "Unrecognized EC fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    radius, atype = FpCommon.parse_fp_token(fp[2..])
    # $stderr << "Radius #{radius} predefined '#{predefined}' ust_atype #{ust_atype}\n"

    radius ||= '3'

    cmd << " -J NCEC#{radius} -R #{radius}"

    atype ||= 'UST:Y'
    cmd << " -P #{atype}"

    cmd
  end
end
