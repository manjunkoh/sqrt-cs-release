classdef NormChanceConstraint
%NORMCHANCECONSTRAINT Container for norm-form chance constraints
%   Represents a chance constraint of the form:
%       P(||x||_2 <= gamma) >= 1 - p
%   where
%       gamma : scalar threshold
%       p     : violation probability in [0,1]
%   The class stores gamma, p and the dimension n (size of x).

    properties
        gamma % scalar threshold
        p     % violation probability (0..1)
        dim   % dimension of x (positive integer)
    end

    methods
        function obj = NormChanceConstraint(gamma, p, n)
            % Constructor: obj = NormChanceConstraint(gamma, p, n)
            if nargin == 0
                obj.gamma = [];
                obj.p = [];
                obj.dim = 0;
                return
            end

            validateattributes(gamma, {'numeric'}, {'scalar'});
            validateattributes(p, {'numeric'}, {'scalar', '>=', 0, '<=', 1});
            validateattributes(n, {'numeric'}, {'scalar', 'integer', 'positive'});

            obj.gamma = double(gamma);
            obj.p = double(p);
            obj.dim = double(n);
        end

        function tf = is_valid(obj)
            tf = ~isempty(obj.gamma) && isscalar(obj.gamma) && isscalar(obj.p) && obj.p >= 0 && obj.p <= 1 && obj.dim > 0;
        end

        function s = toString(obj)
            if isempty(obj.dim) || obj.dim == 0
                s = "NormChanceConstraint (empty)";
            else
                s = sprintf('NormChanceConstraint: dim=%d, gamma=%g, p=%g', obj.dim, obj.gamma, obj.p);
            end
        end

        function [gamma, prob, n] = unpack(obj)
            gamma = obj.gamma;
            prob = obj.p;
            n = obj.dim;
        end

        function constr = toYALMIPConstraint(obj, mu, S)
            %toYALMIPConstraint Build YALMIP constraint for norm chance constraint
            %   constr = mu_norm + sqrt(chi2inv(1-p, n)) * ||S||_2 - gamma <= 0

            if ~obj.is_valid()
                error('NormChanceConstraint:Invalid', 'Constraint object is not valid.');
            end

            prob = obj.p;
            n = obj.dim;

            % Compute chi-square inverse at 1-prob: prefer chi2inv, fallback to gaminv-based formula
            if exist('chi2inv', 'file') == 2
                chi2q = chi2inv(1 - prob, n);
            elseif exist('gaminv', 'file') == 2
                % chi2 with n degrees is Gamma(k=n/2, theta=2)
                chi2q = 2 * gaminv(1 - prob, n/2, 1);
            else
                error('NormChanceConstraint:MissingStats', 'chi2inv or gaminv required from Statistics Toolbox to compute chi-square quantile.');
            end

            % Ensure non-negative
            if chi2q < 0
                chi2q = 0;
            end

            q = sqrt(chi2q);

            % Build expression: norm(mu,2) + q * norm(S,2) - gamma <= 0
            expr = norm(mu, 2) + q * norm(S, 2) - obj.gamma;
            constr = (expr <= 0);
        end
    end
end
