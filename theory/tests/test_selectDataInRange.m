function tests = test_selectDataInRange
% Unit tests for selectDataInRange
    tests = functiontests(localfunctions);
end

function testSelect(testCase)
    t = 0:5;
    I = t.^2;
    [ts, Is] = selectDataInRange(t, I, 2, 4);
    verifyEqual(testCase, ts, 2:4);
    verifyEqual(testCase, Is, [4 9 16]);
end
