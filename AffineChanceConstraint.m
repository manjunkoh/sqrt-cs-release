classdef AffineChanceConstraint
%AFFINECHANCECONSTRAINT Simple container for linear chance constraints
%   Represents a chance constraint of the form:
%       P(a' * x <= b) >= 1 - p
%   where
%       a : vector (alpha)
%       b : scalar (beta)
%       p : violation probability in [0,1]
%       nodes : optional vector of node indices this constraint applies to
%   The class stores a, b, p, the implied size (length of a), and optionally
%   the node indices where this constraint is active.

    properties
        a   % column vector alpha (n x 1)
        b   % scalar beta
        p   % violation probability (scalar in [0,1])
        size % integer n = numel(a)
        nodes % optional: indices of nodes this constraint applies to
    end

    methods
        function obj = AffineChanceConstraint(a, b, p, nodes)
            % Constructor
            %   obj = AffineChanceConstraint(a, b, p)
            %   obj = AffineChanceConstraint(a, b, p, nodes)
            % If no arguments are given, constructs an empty object.
            if nargin == 0
                obj.a = [];
                obj.b = [];
                obj.p = [];
                obj.size = 0;
                obj.nodes = [];
                return
            end

            validateattributes(a, {'numeric'}, {'vector', 'nonempty'});
            validateattributes(b, {'numeric'}, {'scalar'});
            validateattributes(p, {'numeric'}, {'scalar', '>=', 0, '<=', 1});

            obj.a = a(:); % ensure column vector
            obj.b = double(b);
            obj.p = double(p);
            obj.size = numel(obj.a);
            
            if nargin < 4
                obj.nodes = [];
            else
                validateattributes(nodes, {'numeric'}, {'vector', 'integer', 'positive'});
            end
        end

        function tf = is_valid(obj)
            %IS_VALID True if object contains a valid constraint
            tf = ~isempty(obj.a) && isscalar(obj.b) && isscalar(obj.p) && obj.p >= 0 && obj.p <= 1;
        end

        function s = toString(obj)
            %TOSTRING Human-readable description
            if isempty(obj.size) || obj.size == 0
                s = "AffineChanceConstraint (empty)";
            else
                if isempty(obj.nodes)
                    s = sprintf('AffineChanceConstraint: a in R^{%d}, b=%g, p=%g', obj.size, obj.b, obj.p);
                else
                    s = sprintf('AffineChanceConstraint: a in R^{%d}, b=%g, p=%g, nodes=[%s]', ...
                        obj.size, obj.b, obj.p, mat2str(obj.nodes'));
                end
            end
        end

        function [alpha, beta, prob] = unpack(obj)
            %UNPACK Return components
            alpha = obj.a;
            beta = obj.b;
            prob = obj.p;
        end
        
        function node_indices = getNodes(obj)
            %GETNODES Return the node indices for this constraint
            node_indices = obj.nodes;
        end
        
        function [alpha, beta, prob, nodes] = unpackWithNodes(obj)
            %UNPACKWITHNODES Return all components including nodes
            alpha = obj.a;
            beta = obj.b;
            prob = obj.p;
            nodes = obj.nodes;
        end

        function constr = toYALMIPConstraint(obj, mu, S)
            %toYALMIPConstraint Build a YALMIP constraint expression
            %   constr = obj.toYALMIPConstraint(mu, S) returns an inequality
            %   expression suitable for use inside YALMIP when mu and S are
            %   either numeric or YALMIP variables. The returned constraint
            %   represents: mu'*a + Phi^{-1}(1-p) * norm(S'*a) - b <= 0
            %
            %   Requires YALMIP in path if mu or S are sdpvars.

            if ~obj.is_valid()
                error('AffineChanceConstraint:Invalid', 'Constraint object is not valid.');
            end

            % inverse CDF of standard normal at 1-p
            alpha = obj.a;
            beta = obj.b;
            prob = obj.p;

			z = norminv(1 - prob);

            constr = [];
            for k = obj.nodes
                constr = [constr;
                    alpha' * mu(:,k) + z * norm(alpha' * S(:,:,k)) - beta <= 0
                ];
            end
        end
    end
end
