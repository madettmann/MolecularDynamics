classdef Atom < handle
    properties
        id
        pos
        posBefore
        vel
        force
        weight
    end
    methods
        function tf = eq(a,b)
            tf = a.id == b.id;
        end
    end
end