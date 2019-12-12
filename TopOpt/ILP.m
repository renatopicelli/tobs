classdef ILP < Optimization

    properties
        epsilons;
        Values;
        RelaxedLimits;  % Upper Bounds on Constraints after relaxation
        VariableLimits; % Bounds on Design Variables
        FlipLimits; % Limits vector for flipping the variables

    end

    methods

        function self = ILP (epsilons, Limits, Values, DesignVariables, FlipLimits, Optimize)
            if (nargin < 6)
                Optimize = 'Minimize';
            end
            self = self@Optimization (Limits, DesignVariables, Optimize);
            self.epsilons = epsilons;
            self.Values = Values;
            self.FlipLimits = FlipLimits;
%             self.RelaxedLimits = RelaxConstraints (self);
            self.VariableLimits = GetVariableLimits (self);
        end

        function Limits = RelaxConstraints (self)
            Targets = self.Limits - self.Values;
            Limits = Targets;
            for ind = 1:self.nConstraints
                Fraction = self.epsilons(ind) * abs(self.Values(ind));
                if (Targets(ind) > Fraction)
                    Limits(ind) = Fraction;
                elseif (Targets(ind) < -Fraction)
                    Limits(ind) = -Fraction;
                end
            end
        end

        function [Limits, Sensitivities, Obj] = NormalizationRelax (self, Sensitivities, Obj)
            nD = size(Sensitivities, 1);
            Norms = max(max(abs(Sensitivities)), eps);
            d = max(abs(Obj));
            Obj = Obj/max(d,eps);
            Sensitivities = Sensitivities ./ repmat (Norms, nD, 1);
            NormedTargets = (self.Limits - self.Values) ./ Norms';
            AllowedChanges = self.epsilons .* abs(self.Values) ./ Norms';
            A = (NormedTargets > AllowedChanges);
            B = (NormedTargets < -AllowedChanges);
            Limits = (A - B) .* AllowedChanges + (1 - (A + B)) .* NormedTargets;
%             Limits = (A + B) .* AllowedChanges + (1 - (A + B)) .* NormedTargets;
%             Sensitivities(:,B) = -Sensitivities(:,B);
        end

        function VarLimits = GetVariableLimits (self)
            VarLimits = zeros(self.nDesignVariables, 2);
            VarLimits(:, 1) = -1.0 * (abs(self.DesignVariables - 1) < 0.001);
            VarLimits(:, 2) = 1.0 * (abs(self.DesignVariables) < 0.001);
        end

        function [UpdatedVariables, PythonObjCoeff, ...
		    PythonConstCoeff, PythonRelaxedLimits, ...
		    PythonLowerLimits, PythonUpperLimits, PythonnDesignVariables] = Optimize (self, ObjCoeff, ConstCoeff, options)
            [self.RelaxedLimits, ConstCoeff, ObjCoeff] = NormalizationRelax (self, ConstCoeff, ObjCoeff);
            exitflag = 0;
            if (~self.isMinimize)
                ObjCoeff = -ObjCoeff;
            end
            if (nargin < 4) %number of inputs (arguments)
                options.Symmetry = 0;
                options.Optimizer = 'cplex';
            else
               if (~isfield(options, 'Symmetry'))
                   options.Symmetry = 0;
               end
               if (~isfield(options, 'Optimizer'))
                   options.Optimizer = 'cplex';
               end
            end
            Symmetry = options.Symmetry;
            Optimizer = options.Optimizer;
            DesignVariables = self.DesignVariables;
            if (self.SpaceDimension == 2)
                if (Symmetry == 1)
                    C = zeros(self.nDesignVariables / 2, self.nConstraints);
                    Xend = self.Elements(1);
                    Yend = self.Elements(2) / 2;
                    ObjCoeff = reshape(ObjCoeff, self.Elements);
                    ObjCoeff(1:Xend, 1:Yend) = ObjCoeff(1:Xend, 1:Yend) + fliplr(ObjCoeff(:, Yend+1:self.Elements(2)));
                    for ind = 1:self.nConstraints
                        CC = ConstCoeff(:, ind);
                        CC = reshape(CC, self.Elements);
                        CC = CC(1:Xend, 1:Yend) + fliplr(CC(:, Yend+1:self.Elements(2)));
                        C(:, ind) = CC(:);
                    end
                    ConstCoeff = C;
                elseif (Symmetry == 2)
                    C = zeros(self.nDesignVariables / 4, self.nConstraints);
                    Xend = self.Elements(1) / 2;
                    Yend = self.Elements(2) / 2;
                    ObjCoeff = reshape(ObjCoeff, self.Elements);
                    ObjCoeff(1:Xend, 1:Yend) = ObjCoeff(1:Xend, 1:Yend) + fliplr(ObjCoeff(1:Xend, Yend+1:self.Elements(2))) + ...
                                               flipud(ObjCoeff(Xend+1:self.Elements(1), 1:Yend)) + ...
                                               rot90(ObjCoeff(Xend+1:self.Elements(1), Yend+1:self.Elements(2)), 2);
                    for ind = 1:self.nConstraints
                        CC = ConstCoeff(:, ind);
                        CC = reshape(CC, self.Elements);
                        CC = CC(1:Xend, 1:Yend) + fliplr(CC(1:Xend, Yend+1:self.Elements(2))) + ...
                                               flipud(CC(Xend+1:self.Elements(1), 1:Yend)) + ...
                                               rot90(CC(Xend+1:self.Elements(1), Yend+1:self.Elements(2)), 2);
                        C(:, ind) = CC(:);
                    end
                    ConstCoeff = C;
                else
                    Xend = self.Elements(1);
                    Yend = self.Elements(2);
                    ObjCoeff = reshape(ObjCoeff, self.Elements);
                end
            else % SpaceDimension = 3
                if (Symmetry == 3)
                    C = zeros(self.nDesignVariables/8, self.nConstraints);
                    Xend = self.Elements(1) / 2;
                    Yend = self.Elements(2) / 2;
                    Zend = self.Elements(3) / 2;
                    ObjCoeff = reshape(ObjCoeff, self.Elements);
                    ObjCoeff(1:Xend, 1:Yend, 1:Zend) = ObjCoeff(1:Xend, 1:Yend, 1:Zend) + ...
                                                       fliplr(ObjCoeff(1:Xend, Yend+1:self.Elements(2), 1:Zend)) + ...
                                                       flipud(ObjCoeff(Xend+1:self.Elements(1), 1:Yend, 1:Zend)) + ...
                                                       rot90(ObjCoeff(Xend+1:self.Elements(1), Yend+1:self.Elements(2), 1:Zend), 2) + ...
                                                       flip(ObjCoeff(1:Xend, 1:Yend, Zend+1:self.Elements(3)), 3) + ...
                                                       flip(fliplr(ObjCoeff(1:Xend, Yend+1:self.Elements(2), Zend+1:self.Elements(3))), 3) + ...
                                                       flip(flipud(ObjCoeff(Xend+1:self.Elements(1), 1:Yend, Zend+1:self.Elements(3))), 3) + ...
                                                       flip(rot90(ObjCoeff(Xend+1:self.Elements(1), Yend+1:self.Elements(2), Zend+1:self.Elements(3)), 2), 3);
                    for ind = 1:self.nConstraints
                        CC = ConstCoeff(:, ind);
                        CC = reshape(CC, self.Elements);
                        CC = CC(1:Xend, 1:Yend, 1:Zend) + ...
                             fliplr(CC(1:Xend, Yend+1:self.Elements(2), 1:Zend)) + ...
                             flipud(CC(Xend+1:self.Elements(1), 1:Yend, 1:Zend)) + ...
                             rot90(CC(Xend+1:self.Elements(1), Yend+1:self.Elements(2), 1:Zend), 2) + ...
                             flip(CC(1:Xend, 1:Yend, Zend+1:self.Elements(3)), 3) + ...
                             flip(fliplr(CC(1:Xend, Yend+1:self.Elements(2), Zend+1:self.Elements(3))), 3) + ...
                             flip(flipud(CC(Xend+1:self.Elements(1), 1:Yend, Zend+1:self.Elements(3))), 3) + ...
                             flip(rot90(CC(Xend+1:self.Elements(1), Yend+1:self.Elements(2), Zend+1:self.Elements(3)), 2), 3);
                        C(:, ind) = CC(:);
                    end
                    ConstCoeff = C;
                else
                    Xend = self.Elements(1);
                    Yend = self.Elements(2);
                    Zend = self.Elements(3);
                    ObjCoeff = reshape(ObjCoeff, self.Elements);
                end
            end
            DesignVariables = reshape(DesignVariables, self.Elements);
            LowerLimits = reshape(self.VariableLimits(:, 1), self.Elements);
            UpperLimits = reshape(self.VariableLimits(:, 2), self.Elements);
            if (self.SpaceDimension == 2)
                DesignVariables = DesignVariables(1:Xend, 1:Yend);
                DesignVariables = DesignVariables(:);
                ObjCoeff = ObjCoeff(1:Xend, 1:Yend);
                ObjCoeff = ObjCoeff(:);
                LowerLimits = LowerLimits(1:Xend, 1:Yend);
                LowerLimits = LowerLimits(:);
                UpperLimits = UpperLimits(1:Xend, 1:Yend);
                UpperLimits = UpperLimits(:);
            else
                DesignVariables = DesignVariables(1:Xend, 1:Yend, 1:Zend);
                DesignVariables = DesignVariables(:);
                ObjCoeff = ObjCoeff(1:Xend, 1:Yend, 1:Zend);
                ObjCoeff = ObjCoeff(:);
                LowerLimits = LowerLimits(1:Xend, 1:Yend, 1:Zend);
                LowerLimits = LowerLimits(:);
                UpperLimits = UpperLimits(1:Xend, 1:Yend, 1:Zend);
                UpperLimits = UpperLimits(:);
            end
            nDesignVariables = length(DesignVariables);
            C1 = (DesignVariables == 0);
            C2 = -(DesignVariables == 1);
            Trun = C1 + C2;
%             self.RelaxedLimits = [self.RelaxedLimits; self.FlipLimits(1) * nDesignVariables; ...
%                                     self.FlipLimits(2) * nDesignVariables];
            self.RelaxedLimits = [self.RelaxedLimits; self.FlipLimits * nDesignVariables];
%             ConstCoeff = [ConstCoeff'; C1(:)'; C2(:)'];
            ConstCoeff = [ConstCoeff'; Trun'];
            while (exitflag ~= 1)
                if (strcmpi('intlinprog', Optimizer))
                     OptimizerOptions = optimoptions('intlinprog', 'CutGeneration', 'intermediate', 'RootLPAlgorithm','primal-simplex', ...
                                                     'NodeSelection', 'mininfeas', 'HeuristicsMaxNodes', 100, 'RootLPMaxIter', 60000, ...
                                                     'MaxNodes', 1e5);
%                     ScaleObj = max(abs(ObjCoeff));
%                     ObjCoeff = ObjCoeff/ScaleObj;
%                     ScaleCons = max(abs(ConstCoeff), [], 2);
%                     nc = size(ConstCoeff, 1);
%                     for ind = 1:nc
%                         ConstCoeff(ind, :) = ConstCoeff(ind, :)/ScaleCons(ind);
%                         self.RelaxedLimits(ind) = self.RelaxedLimits(ind)/ScaleCons(ind);
%                     end
                    [x, ObjValue, exitflag, output] = intlinprog (ObjCoeff, 1:nDesignVariables, ...
                                                      ConstCoeff, self.RelaxedLimits, [], [], ...
                                                      LowerLimits, UpperLimits, ...
                                                      OptimizerOptions);
%                 output
		elseif (strcmpi('glpk', Optimizer))
                    % OptimizerOptions = optimoptions('intlinprog', 'CutGeneration', 'intermediate', 'RootLPAlgorithm','primal-simplex', ...
                      %                               'NodeSelection', 'mininfeas', 'HeuristicsMaxNodes', 100, 'RootLPMaxIter', 60000, ...
                      %                               'MaxNodes', 1e5);
%                     ScaleObj = max(abs(ObjCoeff));
%                     ObjCoeff = ObjCoeff/ScaleObj;
%                     ScaleCons = max(abs(ConstCoeff), [], 2);
%                     nc = size(ConstCoeff, 1);
%                     for ind = 1:nc
%                         ConstCoeff(ind, :) = ConstCoeff(ind, :)/ScaleCons(ind);
%                         self.RelaxedLimits(ind) = self.RelaxedLimits(ind)/ScaleCons(ind);
%                     end
		    param.msglev = 1;
		    param.btrack = 1;
		    % param.dual = 1;
		    % param.itlim = 100;
		    % param.branch = 5;
		    % param.btrack = 1

		    % param.tmlim = 2000; %ms
                    [x, ObjValue, exitflag, output] = glpk(ObjCoeff, ConstCoeff, self.RelaxedLimits, LowerLimits, UpperLimits, "UU", repmat("I",[1,nDesignVariables]), 1, param );
		    exitflag = 1

		elseif (strcmpi('python_cplex', Optimizer))
                    PythonObjCoeff = ObjCoeff;
		    PythonConstCoeff = ConstCoeff;
		    PythonRelaxedLimits = self.RelaxedLimits;
		    PythonLowerLimits = LowerLimits;
		    PythonUpperLimits = UpperLimits;
		    PythonnDesignVariables = nDesignVariables;
		    ObjValue = 1;
		    exitflag = 1;
		    x = ones(nDesignVariables,1);

                else
                    OptimizerOptions = cplexoptimset('cplex');
                    OptimizerOptions.mip.strategy.variableselect = 3; % OptimizerOptions.BranchStrategy = 'strong';
                    OptimizerOptions.mip.strategy.nodeselect = 2; %OptimizerOptions.Display = 'on';
                    [x, ObjValue, exitflag, output] = cplexmilp(ObjCoeff, ConstCoeff, ...
                                                      self.RelaxedLimits, [], [], [], [], [], ...
                                                      LowerLimits, UpperLimits, ...
                                                      repmat('I', 1, nDesignVariables), [], OptimizerOptions);
%                 output
                end
                self.RelaxedLimits = self.RelaxedLimits * 0.3;
                if (ObjValue == 0)
                    exitflag = 0;
                    self.RelaxedLimits = self.RelaxedLimits * 4;
                    disp('Struck ...');
		    exit;
                end
            end
            disp (['Optimized Objective Value from ILP = ' num2str(ObjValue, 10)]);
            x = round (x);
            UpdatedVariables = DesignVariables + x;
            if (self.SpaceDimension == 2)
                UpdatedVariables = reshape(UpdatedVariables, Xend, Yend);
                if (Symmetry == 1)
                    UpdatedVariables = [UpdatedVariables fliplr(UpdatedVariables)];
                elseif (Symmetry == 2)
                    UpdatedVariables = [UpdatedVariables fliplr(UpdatedVariables)];
                    UpdatedVariables = [UpdatedVariables ; flipud(UpdatedVariables)];
                end
            else % SpaceDimension = 3
                UpdatedVariables = reshape(UpdatedVariables, Xend, Yend, Zend);
                if (Symmetry == 3)
                    UpdatedVariables = [UpdatedVariables fliplr(UpdatedVariables)];
                    UpdatedVariables = [UpdatedVariables; flipud(UpdatedVariables)];
                    UpdatedVariables(:, :, Zend+1:self.Elements(3)) = flip(UpdatedVariables, 3);
                end
            end
%             UpdatedVariables = round (self.DesignVariables + x);
        end

    end % end mehtods

end
