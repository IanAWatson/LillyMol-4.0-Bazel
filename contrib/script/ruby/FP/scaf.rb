# frozen_string_literal: true

# SCAF class for gfp_make

require_relative 'lib/fp_common'

# Scaffold class.
class SCAF
  attr_reader :description

  def initialize
    @rx = Regexp.new('^SCAF')
    @description = 'Scaffold fingerprints'
    @executable = 'molecular_abstractions'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^SCAF/.match(fp)
    raise "Unrecognized SCAF fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)

    if fp.include?('I')
      cmd << " -a 'scaffold(ISO FP)'"
    else
      cmd << " -a 'scaffold(FP)'"
    end
    cmd
  end
end
