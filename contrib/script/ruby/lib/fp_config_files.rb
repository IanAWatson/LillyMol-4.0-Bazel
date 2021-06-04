# frozen_string_literal: true

# Part of gfp_make
# Return a map of fingerprint names to object that can process
# that kind of fingerprint.
# By convention, the name of the class is the uppercase of
# the config file name.

def config_fingerprints(config_dir, verbose)
  fps = {}
  $stderr << "Fetching configs from #{config_dir}\n" if verbose.positive?
  Dir.entries(config_dir).each do |fname|
    next if /^\./.match(fname)
    next if /fp_common.rb/.match(fname)
    next if fname == 'lib'
    $stderr << "Examining #{fname}\n" if verbose > 2
    m = /(\S+)\.rb/.match(fname)
    if not m
      $stderr << "No match file name #{fname}, ignored\n"
      next
    end
    path_name = "#{config_dir}/#{fname}"
    require(path_name)
    class_name = m[1].upcase   # By convention.
    fps[class_name] = eval("#{class_name}.new")
  end

  if verbose.positive?
    $stderr << "#{fps.length} fingerprints recognized\n"
    fps.each do |k, v|
      $stderr << "#{k} #{v} #{v.description}\n"
    end
  end

  fps
end
