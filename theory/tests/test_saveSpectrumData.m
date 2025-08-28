function tests = test_saveSpectrumData
% Unit tests for saveSpectrumData
    tests = functiontests(localfunctions);
end

function testFileCreation(testCase)
    tmp = tempname;
    array1 = [1; 2];
    array2 = [3; 4];
    saveSpectrumData(array1, array2, tmp);
    filename = [tmp, '.txt'];
    verifyTrue(testCase, isfile(filename));
    data = readmatrix(filename);
    verifyEqual(testCase, data, [1 3; 2 4]);
    delete(filename);
end
