function tests = test_saveIFFTPulse
% Unit tests for saveIFFTPulse
    tests = functiontests(localfunctions);
end

function testFileCreation(testCase)
    tmp = tempname;
    array1 = [0; 1];
    array2 = [0.1; 0.2];
    saveIFFTPulse(array1, array2, tmp);
    filename = [tmp, '.txt'];
    verifyTrue(testCase, isfile(filename));
    data = readmatrix(filename);
    verifyEqual(testCase, data, [0 0.1; 1 0.2], 'AbsTol', 1e-12);
    delete(filename);
end
