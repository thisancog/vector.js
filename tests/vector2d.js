'use strict';

let assert = require('assert'),
	Vector2D = require('../dist/vector.min.js');

let vector1 = new Vector2D(0, 0),
	vector2 = new Vector2D(12, 21),
	vector3 = new Vector2D(-5, 9),
	vector4 = new Vector2D(7.4, 1.9),
	vector5 = new Vector2D(-9.2, 0.01);


// the minimum accuracy, e.g. the maximum difference between the computational and the actual result
let epsilon = 0.00001;


/***************************************************************
	Add and addExact
 ***************************************************************/

describe('add() and addExact()', () => {
	it("should add two 2D vectors regardless of order", () => {
		assert.deepEqual(Vector2D.add(vector1, vector2), new Vector2D(12, 21));
		assert.deepEqual(Vector2D.add(vector2, vector1), new Vector2D(12, 21));
		assert.deepEqual(Vector2D.addExact(vector4, vector5), new Vector2D(-1.8, 1.91));
		assert.deepEqual(Vector2D.addExact(vector5, vector4), new Vector2D(-1.8, 1.91));

		// changes vector1!
		vector1.add(vector2);
		assert.deepEqual(vector1, new Vector2D(12, 21));

		// changes vector4!
		vector4.addExact(vector5);
		assert.deepEqual(vector4, new Vector2D(-1.8, 1.91));

		vector1 = new Vector2D(0, 0);
		vector4 = new Vector2D(7.4, 1.9)
	});

	it("should add several 2D vectors", () => {
		assert.deepEqual(Vector2D.add(vector1, vector2, vector3), new Vector2D(7, 30));

		vector2.addExact(vector3, vector4);
		assert.deepEqual(vector2, new Vector2D(14.4, 31.9));
		vector2 = new Vector2D(12, 21);
	});

	it("should return input vector if it is the only one", () => {
		assert.deepEqual(Vector2D.add(vector1), vector1);
	});
});



/***************************************************************
	Subtract and subtractExact
 ***************************************************************/

describe('subtract() and subtractExact()', () => {
	it("should subtract two 2D vectors", () => {
		assert.deepEqual(Vector2D.subtract(vector1, vector2), new Vector2D(-12, -21));
		assert.deepEqual(Vector2D.subtractExact(vector1, vector5), new Vector2D(9.2, -0.01));
	});

	it("should subtract several 2D vectors", () => {
		assert.deepEqual(Vector2D.subtract(vector1, vector2, vector3), new Vector2D(-7, -30));
		assert.deepEqual(Vector2D.subtractExact(vector1, vector4, vector5), new Vector2D(1.8, -1.91));

		vector2.subtract(vector3, vector4);
		assert.deepEqual(vector2, new Vector2D(9.6, 10.1));

		vector4.subtractExact(vector5);
		assert.deepEqual(vector4, new Vector2D(16.6, 1.89));

		vector2 = new Vector2D(12, 21);
		vector4 = new Vector2D(7.4, 1.9);
	});

	it("should return input vector if it is the only one", () => {
		assert.deepEqual(Vector2D.subtract(vector1), vector1);
		assert.deepEqual(Vector2D.subtractExact(vector1), vector1);
	});
});



/***************************************************************
	Scalar multiplication and division
 ***************************************************************/

describe('mult() and div()', () => {
	it("should perform scalar multiplication on a 2D vector", () => {
		assert.deepEqual(Vector2D.mult(vector3, 2), new Vector2D(-10, 18));

		// changes vector3!
		vector3.mult(2);
		assert.deepEqual(vector3, new Vector2D(-10, 18));
	});

	it("should perform scalar division on a 2D vector", () => {
		assert.deepEqual(Vector2D.div(vector3, 2), new Vector2D(-5, 9));

		// changes vector3!
		vector3.div(2);
		assert.deepEqual(vector3, new Vector2D(-5, 9));
	});
});



/***************************************************************
	Length, normalise and limit
 ***************************************************************/

describe('length(), normalise() and limit()', () => {
	it("should accurately calculate the length/magnitude of a 2D vector", () => {
		vector3 = new Vector2D(3, -4);
		assert.deepEqual(Vector2D.length(vector3), 5);
	});

	it("should calculate the 2D unit vector", () => {
		let result = Vector2D.normalise(vector2);
		assert.equal(Math.abs(result.x - 0.49613893) < epsilon, true);
		assert.equal(Math.abs(result.y - 0.86824314) < epsilon, true);
	});
});
