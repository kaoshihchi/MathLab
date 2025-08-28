function tests = test_GaussianAmp
% Unit tests for GaussianAmp
    tests = functiontests(localfunctions);
end

function testValue(testCase)
    y = GaussianAmp(0, 1, 0, 1, 0);
    verifyEqual(testCase, y, 1, 'AbsTol', 1e-12);
end
