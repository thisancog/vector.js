'use strict';

let assert = require('assert'),
	Vector2D = require('../dist/vector.min.js');

let vector1 = new Vector2D(0, 0),
	vector2 = new Vector2D(12, 21),
	vector3 = new Vector2D(-5, 9),
	vector4 = new Vector2D(7.4, 1.9),
	vector5 = new Vector2D(-9.2, 0.01);


/***************************************************************
	Add and addExact
 ***************************************************************/

describe('Vector2D.add() and Vector2D.addExact()', () => {
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
	});

	it("should add several 2D vectors", () => {
		assert.deepEqual(Vector2D.add(vector1, vector2, vector3), new Vector2D(19, 51));

		// changes vector2!
		vector2.addExact(vector3, vector4);
		assert.deepEqual(vector2, new Vector2D(5.2, 31.91));
	});

	it("should return input vector if it is the only one", () => {
		assert.deepEqual(Vector2D.add(vector1), vector1);
	});
});



/***************************************************************
	Subtract and subtractExact
 ***************************************************************/

describe('Vector2D.subtract() and Vector2D.subtractExact()', () => {
	it("should subtract two 2D vectors", () => {
		assert.deepEqual(Vector2D.subtract(vector1, vector2), new Vector2D(6.8, -10.91));
		assert.deepEqual(Vector2D.subtractExact(vector1, vector5), new Vector2D(21.2, 20.99));
	});

	it("should subtract several 2D vectors", () => {
		assert.deepEqual(Vector2D.subtract(vector1, vector2, vector3), new Vector2D(11.8, -19.91));
		assert.deepEqual(Vector2D.subtractExact(vector1, vector4, vector5), new Vector2D(23, 19.08));

		// changes vector2!
		vector2.subtract(vector3, vector4);
		assert.deepEqual(vector2, new Vector2D(12, 21));

		// changes vector4!
		vector4.subtractExact(vector5);
		assert.deepEqual(vector4, new Vector2D(7.4, 1.9));
	});

	it("should return input vector if it is the only one", () => {
		assert.deepEqual(Vector2D.subtract(vector1), vector1);
		assert.deepEqual(Vector2D.subtractExact(vector1), vector1);
	});
});


