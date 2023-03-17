
function beta(x)
    _beta = (9 / 4) - (3 / 2) * (x^ 2) + ((22 / 3)) * x - (7 / 3) * x^ 3
    return _beta
end



function _gamma(x)
    chunk = ((161 / 36)) * (1 / 2) * x^ 2
    chunk += (-(109 / 24)) * (1 / 3) * x^ 4
    chunk += (5 / 18) * x^ 6
    return chunk / 4
end


function pthree(b, d, a)
    return (-b +sqrt(d)) / (2 * a)
end


function SixPointDelta(x)
    alpha = 28
    if x > -3  && x <= -2
        x = x + 3
        tempb = beta(x)
        tempgamma = _gamma(x)
        d =tempb^2 - 4 * alpha * tempgamma
        pm3 = pthree(tempb, d, alpha)
        return pm3
    elseif x > -2 && x <= -1
        x = x + 2
        tempb = beta(x)
        tempgamma = _gamma(x)
        d = tempb^ 2 - 4 * alpha * tempgamma
        pm3 = pthree(tempb, d, alpha)
        pm2 = -3 * pm3 - (1 / 16) + (1 / 8) * (x^2) + (1 / 12) * (- 1) * x + (1 / 12) * x^2
        return pm2
    elseif x > -1 && x <= 0
        x = x + 1
        tempb = beta(x)
        tempgamma = _gamma(x)
        d =tempb^2 - 4 * alpha * tempgamma
        pm3 = pthree(tempb, d, alpha)
        pm1 = 2 * pm3 + (1 / 4) + (1 / 6) * (4) * x - (1 / 6) * x^3
        return pm1

    elseif x > 0 && x <= 1
        tempb = beta(x)
        tempgamma = _gamma(x)
        d = tempb^2 - 4 * alpha * tempgamma
        pm3 = pthree(tempb, d, alpha)
        p = 2 * pm3 + (5 / 8) - (1 / 4) * (x^2)
        return p
    elseif x > 1 && x <= 2
        x = x - 1
        tempb = beta(x)
        tempgamma = _gamma(x)
        d = tempb^ 2 - 4 * alpha * tempgamma
        pm3 = pthree(tempb, d, alpha)
        pp1 = -3 * pm3 + (1 / 4) - (1 / 6) * (4) * x + (1 / 6) * x^ 3
        return pp1
    elseif x > 2 && x <= 3
        x = x - 2
        tempb = beta(x)
        tempgamma = _gamma(x)
        d = tempb^ 2 - 4 * alpha * tempgamma
        pm3 = pthree(tempb, d, alpha)
        pp2 = pm3 - (1 / 16) + (1 / 8) * (x^ 2) - (1 / 12) * (- 1) * x - (1 / 12) * x^ 3
        return pp2
    end
    return 0
    
end



