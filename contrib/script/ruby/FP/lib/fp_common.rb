# Functions common across different fingerprint types.

module FpCommon

extend self

def initial_command_stem(executable, first_in_pipeline:, extra_qualifiers:)
  cmd = "#{executable}"
  if first_in_pipeline
    cmd << " -A 2 -g all -l"
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
def parse_fp_token_inner(fp)
  return nil, nil if fp.empty?
  $stderr << "parse_fp_token_inner processing #{fp}\n"
  f = fp.split(':')
  ust_atype = if f.length == 1
               nil
             else
               f[1]
             end
  m = /^([A-Z]*)(\d*)/.match(f[0])
  predefined = if m[1].empty?
                 nil
               else
                 m[1]
               end

  radius = if m[2].empty?
             nil
           else
              m[2]
           end

  if predefined && ust_atype
    raise "Cannot specify both predefined #{predefined} and ust_atype #{ust_atype} from #{fp}"
  end

  return radius, predefined if predefined
  return radius, "UST:#{ust_atype}" if ust_atype
  return radius, nil
end

# Detect whether `fp` contains the token ':fixed'. If so, remove it and pass what remains
# to parse_fp_token_inner.
# Returns whatever parse_fp_token_inner returns, together with a boolean for presence of :fixed

def parse_fp_token(fp)
  return nil, nil, false if fp.empty?

  if fp.match(/:fixed/)
    fp.slice!(':fixed')
    return *parse_fp_token_inner(fp), true
  elsif fp.match(/:sparse/)
    return *parse_fp_token_inner(fp), false
  else
    return *parse_fp_token_inner(fp), false
  end
end

end
