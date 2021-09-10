# frozen_string_literal: true

# EC fingerprint module for gfp_make

require_relative 'lib/fp_common'

# EC fingerprint.
class NEC
  attr_reader :description

  def initialize
    @rx = Regexp.new('^NEC')
    @description = 'Circular fingerprints (new)'
    @executable = 'ec_fingerprint'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^NEC([A-Z]*)(\d+)*:*(\S+)*/.match(fp)
    raise "Unrecognized NEC fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    radius, atype = FpCommon.parse_fp_token(fp[3..])
    # $stderr << "Radius #{radius} predefined '#{predefined}' ust_atype #{ust_atype}\n"

    radius ||= '3'

    cmd << " -J NCEX#{radius}" unless /-J /.match(extra_qualifiers)
    cmd << " -R #{radius}" unless /-R \d/.match(extra_qualifiers)

    atype ||= 'UST:Y'
    cmd << " -P #{atype}"

    cmd
  end
end
