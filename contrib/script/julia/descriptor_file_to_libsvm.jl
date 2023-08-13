# Converts an activity file and a descriptor file to libsvm

using ArgMacros
using DelimitedFiles

function read_activity_file(fname::String)::Dict{String, String}
  result = Dict{String, String}()
  raw_data = DelimitedFiles.readdlm(fname, ' ', String, skipstart=1)
  for row in eachrow(raw_data)
    result[row[1]] = row[2]
  end
  return result
end

function main()
  @inlinearguments begin
    @argumentrequired String activity_file "-A" "--activity_file"
    @argumentflag verbose "-v" "--verbose"
    @positionalrequired String csv_file
  end

  activity = read_activity_file(activity_file)
  if verbose
    print(stderr, "Activity file $(activity_file) contains $(length(activity)) rows\n")
  end

  descriptors = DelimitedFiles.readdlm(csv_file, ' ', String, skipstart=1)
  for row in eachrow(descriptors)
    id = row[1]
    act = get(activity, id, "")
    if act == ""
      print(stderr, "No activity for $(id)\n")
      exit(1)
    end
    print(act)
    for i in 2:length(row)
      print(" $(i-2):$(row[i])")
    end
    print("\n")
  end
end

main()
