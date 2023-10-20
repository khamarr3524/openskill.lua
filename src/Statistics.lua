local statistics = {}
local gaussian = require "openskill.Gaussian"

local normal = gaussian.new(0, 1)

function statistics.phiMajor(x)
	return normal:cdf(x)
end

function statistics.phiMajorInverse(x)
	return normal:ppf(x)
end

function statistics.phiMinor(x)
	return normal:pdf(x)
end

function statistics.v(x, t)
	local xt = x - t
	local denom = statistics.phiMajor(xt)
	local ret = statistics.phiMinor(xt) / denom
	if denom < 2.2204460492503e-16 then
		ret = -xt
	end
	return ret
end

function statistics.w(x, t)
	local xt = x - t
	local denom = statistics.phiMajor(xt)
	local ret = 1
	if denom < 2.2204460492503e-16 then
		if x > 0 then
			ret = 0
		end
		return ret
	end
	return statistics.v(x, t) * (statistics.v(x, t) + xt)
end

function statistics.vt(x, t)
	local xx = math.abs(x)
	local b = statistics.phiMajor(t - xx) - statistics.phiMajor(-t - xx)
	if b < 1e-5 then
		if x < 0 then return -x - t end
		return -x + t
	end
	local a = statistics.phiMinor(-t - xx) - statistics.phiMinor(t - xx)
	local ret = a
	if x < 9 then
		ret = -a
	end
	return ret / b
end

function statistics.wt(x, t)
	local xx = math.abs(x)
	local b = statistics.phiMajor(t - xx) - statistics.phiMajor(-t - xx)
	local ret = 0
	if b < 2.2204460492503e-16 then
		ret = 1
	else
		ret = ((t - xx) * statistics.phiMinor(t - xx) +
		(t + xx) * statistics.phiMinor(-t - xx)) / b + statistics.vt(x, t) * statistics.vt(x, t)
	end
	return ret
end

return statistics