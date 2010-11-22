function add_task(handle, vertex, functionname, priority) %!eml
    eml.ceval('emx_add_task', handle, vertex, functionname, priority);
end
