# frozen_string_literal: true

# W class for gfp_make

require_relative 'lib/fp_common'

# Class for W topological torsion fingerprints.
class W
  attr_reader :description

  def initialize
    @rx = Regexp.new('^W')
    @description = 'iwdescr fingerprints'
    @executable = 'iwdescr'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^W(\d+)*/.match(fp)
    raise "Unrecognized TT fp form '#{fp}'" unless m

    cmd = @executable.dup
    cmd << " -G FILTER" unless first_in_pipeline
    cmd << " -G #{m[1]}" if m[1].to_i > 1
    cmd << " -G ALL"
    # path_length, atype, fixed = FpCommon.parse_fp_token(fp[2..])

    cmd
  end
end
