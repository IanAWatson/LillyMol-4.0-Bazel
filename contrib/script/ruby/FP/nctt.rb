# frozen_string_literal: true

# TT class for gfp_make

require_relative 'lib/fp_common'

# Class for TT topological torsion fingerprints.
class NCTT
  attr_reader :description

  def initialize
    @rx = Regexp.new('^NCTT')
    @description = 'Topological torsion fingerprint'
    @executable = 'topological_torsion'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^NCTT(\d+)*/.match(fp)
    raise "Unrecognized TT fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    path_length, atype = FpCommon.parse_fp_token(fp[2..])

    cmd << ' -J NCTT'
    cmd << " -P #{atype}" if atype
    cmd
  end
end
