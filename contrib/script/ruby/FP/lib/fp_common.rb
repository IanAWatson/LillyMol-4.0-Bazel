# Functions common across different fingerprint types.

module FpCommon

extend self

def initial_command_stem(executable, first_in_pipeline:, extra_qualifiers:)
  cmd = "#{executable}"
  if first_in_pipeline
    cmd << " -g all -l"
  else
    cmd << " -f"
  end

  if extra_qualifiers
    cmd << " #{extra_qualifiers}"
  end

  cmd
end

def legacy_atom_type_specification(legacy_atom_type)
  case legacy_atom_type
  when 'C'
    return '-P complex'
  when 'SB'
    return '-P sb'
  when 'TT'
    return '-P tt'
  else
    raise "Unrecognized legacy atom type #{legacy_atype}"
  end
end

# Parse something that looks like
#  '' (empty)
#  3
#  C
#  C3
#  C:xxx
#  C3:xxx
#  3:xxx
#  :xxx
#  Return a tuple of
#    a. the 3 (the radius)
#    b. the atom type. If derived from xxx, then UST:xxx
#  Raise an exception if both atom type forms are specified.
def parse_fp_token(fp)
  return nil, nil if fp.empty?
  f = fp.split(':')
  if f.length == 1
    ust_atype = nil
  else
    ust_atype = f[1]
  end
  m = /^([A-Z]*)(\d*)/.match(f[0])
  if m[1].empty?
    predefined = nil
  else
    predefined = m[1]
  end
  if m[2].empty?
    radius = nil
  else
    radius = m[2]
  end

  if predefined && ust_atype
    raise "Cannot specify both predefined #{predefined} and ust_atype #{ust_atype} from #{fp}"
  end

  return radius, predefined if predefined
  return radius, "UST:#{ust_atype}" if ust_atype
  return radius, nil
end

end
