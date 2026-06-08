module Fileio
"""
Misc file import/export operations
"""

    export sort_steps


    function sort_steps(filename)
        """
        Sort steps from alphabetical/lexicographic to natural
        """
        if startswith(filename, "step_")
            return parse(Int, split(split(filename, "_")[2], "_")[1])
        else
            return -1  # So that non-step files are sorted to the start
        end
    end

end