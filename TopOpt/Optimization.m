classdef Optimization

    properties
        isMinimize;        
        DesignVariables;
        nDesignVariables;
        Elements;
        SpaceDimension;
        nConstraints;
        Limits; % Assuming that all constraints are $\leq$ type. Any other type constraints 
                % need to be converted to this type.
    end
    
    methods
        
        function self = Optimization (Limits, DesignVariables, Optimize)
            if (nargin < 2)
                self.isMinimize = true;
            else
                self.isMinimize = strcmpi (Optimize, 'Minimize');
            end
            self.Elements = size(DesignVariables);
            self.SpaceDimension = length(self.Elements); % HOW ABOUT MACRO?
            self.DesignVariables = DesignVariables(:);
            self.nDesignVariables = prod(self.Elements);
            self.Limits = Limits;
            self.nConstraints = length (Limits);
        end
        
    end
end