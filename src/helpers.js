var isNumeric = function(...numbers) {
	return numbers.filter(n => typeof n === 'undefined' || Array.isArray(n) || isNaN(parseFloat(n)) || !isFinite(n)).length === 0;
}

var sumExact = function(...summands) {
	if (summands.length == 0) return;
	if (summands.length == 1) return summands[0];

	let sum = summands[0],
		compensation = 0.0,
		decimals = 0;

	for (let i = 1; i < summands.length; i++) {
		let summand = summands[i],
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