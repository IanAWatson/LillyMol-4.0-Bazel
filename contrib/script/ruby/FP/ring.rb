# frozen_string_literal: true

# RING class for gfp_make

require_relative 'lib/fp_common'

# Ring class.
class RING
  attr_reader :description

  def initialize
    @rx = Regexp.new('^RING')
    @description = 'Ring atom fingerprints'
    @executable = 'molecular_abstractions'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^RING/.match(fp)
    raise "Unrecognized RING fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)

    cmd << " -c -a 'rings(FP)'"
    cmd
  end
end
