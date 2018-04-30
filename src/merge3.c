// arr1: input=inf end-times
// arr4: output=idx of inf end-times in the merged times
// arr2: input=inf start-times, and/or bolus time & sampl
// arr3: input=idx of inf start-times, and/or bolus time before merge;
//      output=idx inf start-times, and/or bolus time after merge
// arr5: output=overallmerged times
void mergeArrays(double *arr1, double *arr2, int *arr3, int *arr4, double *arr5, int *n1, int *n2, int *n3)
{
    int i = 0, j = 0, k = 0, l = 0, m = 0;

    // Traverse both array
    while (i<*n1 && j <*n2)
    {
        // Check if current element of first
        // array is smaller than current element
        // of second array. If yes, store first
        // array element and increment first array
        // index. Otherwise do same with second array
        if (arr1[i] <= arr2[j]) {
            arr5[l] = arr1[i];
            arr4[i++] = ++l;
        } else {
        // Check if a dosing time
            if (k<*n3 && arr3[k] == m) arr3[k++] = l+1;
            arr5[l++] = arr2[j++];
            m++;
		}
    }

    // Store remaining elements of first array
    while (i < *n1) {
        arr5[l] = arr1[i];
        arr4[i++] = ++l;
    }
    // Store remaining elements of second array
    while (j < *n2) {
        // Check if a dosing time
        if (k<*n3 && arr3[k] == m) arr3[k++] = l+1;
        arr5[l++] = arr2[j++];
        m++;
	}
}

/*
dyn.unload('merge3.dll')
dyn.load('merge3.dll')
arr1 = with(ev[ds,], time + amt/rate)
arr2 = as.double(ev$time)
arr3 = (1:length(ev$evid))[ev$evid>0]-1
s = .C("mergeArrays", arr1, arr2, as.integer(arr3), as.integer(arr1), double(261), length(arr1), length(arr2), length(arr3))
arr1; s[4]; arr3; s[3]; s[5]
arr1-s[[5]][s[[4]]]
s[[2]][arr3+1]-s[[5]][s[[3]]]
*/

