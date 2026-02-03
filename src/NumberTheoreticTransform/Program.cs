/* 
Pi program created by Jason Bock on 2/7/95
This program calculates pi up to a maximum precision of 5000 digits.
A number theoretic transform is used to reduce the complexity of the
multiplication function
*/

/*
The variable M is the modulus number. N is the length of the transform.
limit is used to set the maximum precision for the calculation.
*/

const int M = 50_010_001;
const int halfM = M / 2;
const int N = 5_000;
//short limit;


Console.WriteLine("3.14");
Console.WriteLine(multinv(43_015));
Console.WriteLine(negmod(5_310_431, 415));
Console.WriteLine(negmod(-5_310_431, 415));
Console.WriteLine(mymod(5_310_431, 5_310_431));

/*
The multinv function calculates the multiplicative inverse of ainv using
a variation of the Euclidian algorithm.
*/
static int multinv(int ainv)
{
	int a0, a1, x0, x1, y0, y1, temp_a0, temp_a1, temp_x0, temp_x1, temp_y0, temp_y1;
	int q;

	a0 = M;
	a1 = ainv;
	x0 = 1;
	x1 = 0;
	y0 = 0;
	y1 = 1;

	while (a1 != 0)
	{
		q = a0 / a1;
		temp_a0 = a1;
		temp_a1 = a0 - a1 * q;
		temp_x0 = x1;
		temp_x1 = x0 - x1 * q;
		temp_y0 = y1;
		temp_y1 = y0 - y1 * q;
		a0 = temp_a0;
		a1 = temp_a1;
		x0 = temp_x0;
		x1 = temp_x1;
		y0 = temp_y0;
		y1 = temp_y1;
	}

	return y0;
}

/*
The negmod function calculates x mod y if x is either positive or
negative.
*/

static int negmod(int x, int y)
{
	if (x > 0)
	{
		x = x % y;

		if (x > halfM)
		{
			return x - M;
		}
		else
		{
			return x;
		}
	}
	else
	{
		while (x < halfM)
		{
			x = x + M;
		}

		return x;
	}
}

/*
The mymod function calculates (x * y) mod M. This function is used when
x * y would cause an overflow to occur. Since this function is rather
complex, it is used only when necessary.
*/
static int mymod(int x, int y)
{
	int value, hold;
	long testx, testy, testz;

	/* The additive inverse is used to make both x and y positive. */
	while (x < 0)
	{
		x = x + M;
	}

	while (y < 0)
	{
		y = y + M;
	}

	testx = x;
	testy = y;
	testz = testx * testy;

	if (testz < M)
	{
		value = (int)testz;
	}
	else
	{
		hold = (int)(testz / M);
		value = (int)(testz - M * hold);
	}

	if (value > halfM)
	{
		return (value - M);
	}
	else if (value < -halfM)
	{
		return (value + M);
	}
	else
	{
		return value;
	}
}

/*
The mpmult function takes two vectors, vect1 and vect2 and multiplies
them together using the number theoretic transform.
*/

static void mpmult(int[] vect1, int[] vect2, int[] alp, int[] alpinv, int[] g1)
{
	const int size = 5003;

	int n0, n1, n2, n3, n4, n5, n6, n7, k0, k1, k2, k3, k4, k5, k6, k7, accout, accin, accrt;
	int t, temp, aij;
	int[] x0, x1, h0, h1, g0;

	long check, go = int.MaxValue;
	short i, j;

	x0 = new int[size];
	x1 = new int[size];
	h0 = new int[size];
	h1 = new int[size];
	g0 = new int[size];

	for (i = 0; i < size; i++)
	{
		x1[i] = 0;
		h1[i] = 0;
		x0[i] = 0;
		h0[i] = 0;
		g0[i] = 0;
		g1[i] = 0;
	}

	/*
	The 14 loops of four variables each perform a part of the transform on
	each input vector using the C-T algorithm. When the alp and alpinv arrays
	are accessed, it is determined if the answer would be 1; an addition is then
	perform instead of multiplication. If the element accessed is N/2, then a 
	subtraction occurs. Also, when appropriate a check to determine if an
	overflow will occur is made; this is done by using the go variable. To
	minimize storage area, only two vectors are used and one is cleared out when
	needed.
	*/

	for (k0 = 0; k0 < 5; k0++)
	{
		for (k1 = 0; k1 < 5; k1++)
		{
			for (k2 = 0; k2 < 5; k2++)
			{
				for (k3 = 0; k3 < 5; k3++)
				{
					for (k4 = 0; k4 < 2; k4++)
					{
						for (k5 = 0; k5 < 2; k5++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k6 = 0; k6 < 2; k6++)
								{
									accout = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * k3 + 4 * k4 + 2 * k5 + n0;
									accin = 2500 * k6 + 1250 * k5 + 625 * k4 + 125 * k3 + 25 * k2 + 5 * k1 + k0;
									accrt = (2500 * k6 * n0) % N;

									if (accrt == 0)
									{
										x1[accout] = negmod(vect1[accin] + x1[accout], M);
										h1[accout] = negmod(vect2[accin] + h1[accout], M);
									}
									else if (accrt == 2500)
									{
										x1[accout] = negmod(-vect1[accin] + x1[accout], M);
										h1[accout] = negmod(-vect2[accin] + h1[accout], M);
									}
									else
									{
										// POST-IT
									}
								}
							}
						}
					}
				}
			}
		}
	}
}






