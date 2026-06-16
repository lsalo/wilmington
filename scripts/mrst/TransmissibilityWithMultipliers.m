classdef TransmissibilityWithMultipliers < StateFunction
    % Transmissibility for internal faces. May include an optional
    % pressure-dependent multiplier from a field in the fluid model.
    properties
        
    end
    
    methods
        function pp = TransmissibilityWithMultipliers(model)
            pp@StateFunction(model);
            pp.label = 'T_f';
            assert(isfield(model.operators, 'T'));
            T = value(model.operators.T);
            assert(all(isfinite(T)))
            if any(T < 0)
                warning('Negative transmissibility in %d interfaces', sum(T < 0));
            end
            pp.outputRange = [0, inf];
        end
        
        function T = evaluateOnDomain(sfn, model, state)
            T = model.operators.T.*model.G.faces.tmult(model.operators.internalConn);
        end
    end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
