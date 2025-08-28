function tests = test_superGaussian
% Unit tests for superGaussian
    tests = functiontests(localfunctions);
end

function testValue(testCase)
    y = superGaussian(0, 1, 0, 1, 1, 0);
    verifyEqual(testCase, y, exp(-1), 'AbsTol', 1e-12);
end
