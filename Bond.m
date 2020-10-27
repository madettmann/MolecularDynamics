classdef Bond
    properties
        atom1
        atom2
    end
    methods
        function tf = eq(a,b)
            tf = (a.atom1 == b.atom1 && a.atom2 == b.atom2) || ...
                (a.atom1 == b.atom2 && a.atom2 == a.atom1) || ...
                (a.atom2 == b.atom1 && a.atom1 == b.atom2) || ...
                (a.atom2 == b.atom2 && a.atom1 == b.atom1);
        end
    end
end