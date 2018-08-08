var isNumeric = function(n) {
	return (typeof n !== 'undefined') ? !Array.isArray(n) && !isNaN(parseFloat(n)) && isFinite(n) : undefined;
}

var sumExact = function(data) {
	if (data.length == 0) return;
	if (data.length == 1) return data[0];

	let sum = data[0],
		compensation = 0.0,
		decimals = 0;

	for (let i = 1; i < data.length; i++) {
		let summand = data[i],
			t = sum + summand;

		if ((summand + "").indexOf('.') > -1) {
			let numDec = (summand + "").split(".")[1].length;
			decimals = numDec > decimals ? numDec : decimals;
		}

		compensation += (Math.abs(sum) >= Math.abs(summand)) ? (sum - t) + summand : (summand - t) + sum;
		sum = t;
	}

	let power = Math.pow(10, decimals + 1);
	return Math.round((sum + compensation) * power) / power;
}