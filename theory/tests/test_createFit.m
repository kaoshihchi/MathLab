function tests = test_createFit
% Unit tests for createFit
    tests = functiontests(localfunctions);
end

function testFit(testCase)
    t = linspace(-1e-13, 1e-13, 100)';
    params = [1, 0, 20e-15, 0];
    y = GaussianAmp(t, params(1), params(2), params(3), params(4));
    [fitresult, gof] = createFit(t, y);
    verifyTrue(testCase, isstruct(gof));
    verifyEqual(testCase, fitresult.A, params(1), 'RelTol', 0.5);
end
