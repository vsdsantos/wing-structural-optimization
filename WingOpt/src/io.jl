module InputOutput

using TOML

export read_data, write_data

"""
    write_data(path, data)

Write a string-keyed dictionary of TOML-compatible values to `path`.
"""
function write_data(path::AbstractString, data::AbstractDict{<:AbstractString})
    open(path, "w") do stream
        TOML.print(stream, data)
    end
    return path
end

"""Read a TOML data file and return its contents as a dictionary."""
read_data(path::AbstractString) = TOML.parsefile(path)

end
