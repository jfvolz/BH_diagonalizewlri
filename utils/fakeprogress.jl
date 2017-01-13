# Macro that does nothing. For use when ProgressMeter is not available or not
# wanted.
macro showprogress(x)
    esc(x)
end
