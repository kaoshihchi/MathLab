function tests = test_filterData
% Unit tests for filterData
    tests = functiontests(localfunctions);
end

function testFiltering(testCase)
    I = [1 2 3 4];
    t = [-2 0 1 3];
    [I_f, t_f] = filterData(I, t, 1);
    verifyEqual(testCase, I_f, [1 4]);
    verifyEqual(testCase, t_f, [-2 3]);
end
