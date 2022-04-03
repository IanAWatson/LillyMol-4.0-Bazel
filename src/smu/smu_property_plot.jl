# Pl
using ArgMacros
using CSV
using DataFrames
using Humanize: digitsep
using Plots
using StatsBase

println(STDERR, "Packages loaded")

function plot_distributions(data,    # ::Vector{DataFrame}, 
                            name,    # ::Vector{String},
                            color,   # ::Vector{String},
                            xlabel::String,
                            title::String)
  n = length(data)
  counts = [sum(data[i][:,2]) for i in 1:n]
  normed = [data[i][:,2] / counts[i] for i in 1:n]
  means = []
  for i in 1:n
    w = FrequencyWeights(data[i][:,2])
    push!(means, median(data[i][:,1], w))
  end
  labels = ["$(name[i]) $(digitsep(counts[i])) $(round(means[i], digits=2))" for i in 1:n]
  p = plot(data[1][:,1], normed[1], label=labels[1], color=color[1], legend=:topleft,
        xlabel=xlabel, ylabel="Prevalence", title=title)
  for i in 2:n
    plot!(p, data[i][:,1], normed[i], label=labels[i], color=color[i])
  end
  # display(p)
  savefig(p, "/tmp/smu.png")
end

function main()
  @inlinearguments begin
    @argumentrequired String title "-t" "--title"
    @argumentrequired String xlabel "-l" "--xlabel"
    @argumentrequired String output_stem "-o" "--output_stem"
    @argumentrequired String names_csv "-n" "--name"
    @argumentrequired String colors_csv "-c" "--color"
    @positionalleftover String file_names "files"
  end

  colors = split(colors_csv, ',')
  names = split(names_csv, ',')

  if length(colors) != length(file_names) ||
     length(names) != length(file_names)
    println(STDERR, "Must have same number of colors as files")
    return 1
  end

  data = []
  for fname in file_names
    push!(data, CSV.File(fname, header=1, delim=' ') |> DataFrame)
  end
  plot_distributions(data, names, colors, xlabel, title)
end

main()
