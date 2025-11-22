classdef SCPCode
	%SCPCode Enumeration-like container for SCP status/result codes
	%   Use as: SCPCode.SOLVED, SCPCode.STALLED, SCPCode.REACHED_MAX_ITERS, etc.
	%   Each field is a constant integer value.

	properties (Constant)
		% Success / terminal codes
		SOLVED = int32(0);
		STALLED = int32(1);
		REACHED_MAX_ITERS = int32(2);
		INFEASIBLE = int32(3);
		NUMERICAL = int32(4);
		% Generic error
		OTHER_ERROR = int32(-1);
	end

	methods (Static)
		function name = toString(code)
			%toString Return a short name for a numeric code
			%   name = SCPCode.toString(code)
			switch int32(code)
				case SCPCode.SOLVED
					name = 'SOLVED';
				case SCPCode.STALLED
					name = 'STALLED';
				case SCPCode.REACHED_MAX_ITERS
					name = 'REACHED_MAX_ITERS';
				case SCPCode.INFEASIBLE
					name = 'INFEASIBLE';
				case SCPCode.NUMERICAL
					name = 'NUMERICAL';
				case SCPCode.OTHER_ERROR
					name = 'OTHER_ERROR';
				otherwise
					name = sprintf('UNKNOWN(%d)', int32(code));
			end
		end

		function ok = isValid(code)
			%isValid True if code is one of the defined constants
			val = int32(code);
			ok = any(val == [SCPCode.SOLVED, SCPCode.STALLED, SCPCode.REACHED_MAX_ITERS, SCPCode.INFEASIBLE, SCPCode.NUMERICAL, SCPCode.OTHER_ERROR]);
		end
	end
end