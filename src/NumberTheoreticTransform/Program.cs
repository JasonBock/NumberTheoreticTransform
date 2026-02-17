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

using System.Globalization;
using System.Text;

const int M = 50_010_001;
const int halfM = M / 2;
const int N = 5_000;

Console.WriteLine("3.14");
Console.WriteLine(multinv(43_015));
Console.WriteLine(negmod(5_310_431, 415));
Console.WriteLine(negmod(-5_310_431, 415));
Console.WriteLine(mymod(5_310_431, 5_310_431));

/*
The main part of the program calculates pi using Borwein's quadrati-
cally converging algorithm. For all of the vectors, element 5000 deter-
mines the sign (1 for positive and 0 for negative), element 5001 holds
the exponent of the number and element 5002 holds the size of the
number
*/

int i, aij, q;
var alp = new int[5003];
var alpinv = new int[5003];
var ak = new int[5003];
var bk = new int[5003];
var pk = new int[5003];
var one = new int[5003];
var temp = new int[5003];
var temp2 = new int[5003];
var temp3 = new int[5003];
var ak1 = new int[5003];
var bk1 = new int[5003];
var pk1 = new int[5003];
short correct = 0;
var go = Math.Pow(2, 31) - 1;

/*
Clear out the vectors.
*/

for (i = 0; i < 5003; i++)
{
	alp[i] = 0;
	alpinv[i] = 0;
	ak[i] = 0;
	bk[i] = 0;
	pk[i] = 0;
	one[i] = 0;
	temp[i] = 0;
	temp2[i] = 0;
	temp3[i] = 0;
	ak1[i] = 0;
	bk1[i] = 0;
	pk1[i] = 0;
}

aij = 1;
q = 1;
alp[0] = 1;
alpinv[0] = 1;

/*
This for loop calculates the T and T^-1 elements for the NTT. 
*/

for (i = 1; i < N; i++)
{
	if (Math.Abs(aij) > (go / 12265))
	{
		aij = mymod(aij, 12265);

		while (aij < 0)
		{
			aij = aij + 50010001;
		}
	}
	else
	{
		aij = (aij * 12265) % M;
	}

	if (aij != 1)
	{
		q = multinv(aij);
	}
	else
	{
		q = 1;
	}

	if (aij > 25005000)
	{
		aij = aij - M;
	}
	else if (aij < -25005000)
	{
		aij = aij + M;
	}

	if (q > 25005000)
	{
		q = q - M;
	}
	else if (q < -25005000)
	{
		q = q + M;
	}

	alp[i] = aij;
	alpinv[i] = q;

	while (aij < 0)
	{
		aij = aij + M;
	}

	/*
	ak[], bk[] and pk[] are set to initial values.
	*/

	ak[0] = 2;
	ak[5000] = 1;
	ak[5001] = 0;
	ak[5002] = 1;

	bk[5000] = 1;
	bk[5001] = 0;
	bk[5002] = 1;

	pk[0] = 2;
	pk[5000] = 1;
	pk[5001] = 0;
	pk[5002] = 1;

	one[0] = 1;
	one[5000] = 1;
	one[5001] = 0;
	one[5002] = 1;

	Console.WriteLine("Enter in the limit (max. is 2450)");
	Shared.limit = short.Parse(Console.ReadLine()!, CultureInfo.CurrentCulture);

	/*
	These initial steps set a0 = sqrt(2), b0 = 0 and p0 = 2 + a0.
	*/

	recipsqrt(ak, alp, alpinv, temp);
	Console.WriteLine("Initial Step 1 done");
	mpmult(temp, pk, alp, alpinv, temp2);
	Console.WriteLine("Initial Step 2 done");

	for (i = 0; i < 5003; i++)
	{
		ak[i] = temp2[i];
	}

	Console.WriteLine("Initial Step 3 done");
	mpadd(ak, pk, temp2);
	Console.WriteLine("Initial Step 4 done");

	for (i = 0; i < 5003; i++)
	{
		pk[i] = temp2[i];
	}

	Console.WriteLine("Initial Step 5 done");

	/*
	This while loop performs the pi iteration
	*/

	while (correct < Shared.limit)
	{
		recipsqrt(ak, alp, alpinv, temp);
		Console.WriteLine($"Step A done at correct = {correct}");

		for (i = 0; i < 5003; i++)
		{
			temp2[i] = temp[i];
		}

		Console.WriteLine($"Step B done at correct = {correct}");
		mpmult(ak, temp2, alp, alpinv, temp3);
		Console.WriteLine($"Step C done at correct = {correct}");
		mpadd(temp, temp3, ak1);
		Console.WriteLine($"Step D done at correct = {correct}");
		ptfive(ak1, temp2);
		Console.WriteLine($"Step E done at correct = {correct}");

		for (i = 0; i < 5003; i++)
		{
			ak1[i] = temp2[i];
		}

		Console.WriteLine($"Step F done at correct = {correct}");
		mpadd(one, bk, temp);
		Console.WriteLine($"Step G done at correct = {correct}");
		mpmult(temp, temp3, alp, alpinv, bk1);

		for (i = 0; i < 5003; i++)
		{
			temp3[i] = bk1[i];
		}

		Console.WriteLine($"Step H done at correct = {correct}");
		mpadd(ak, bk, temp);
		Console.WriteLine($"Step I done at correct = {correct}");
		recip(temp, alp, alpinv, temp2);
		Console.WriteLine($"Step J done at correct = {correct}");
		mpmult(temp3, temp2, alp, alpinv, bk1);
		Console.WriteLine($"Step K done at correct = {correct}");
		mpadd(one, ak1, temp);
		Console.WriteLine($"Step L done at correct = {correct}");
		mpmult(bk1, temp, alp, alpinv, temp2);
		Console.WriteLine($"Step M done at correct = {correct}");
		mpmult(pk, temp2, alp, alpinv, temp3);
		Console.WriteLine($"Step N done at correct = {correct}");
		mpadd(one, bk1, temp);
		Console.WriteLine($"Step O done at correct = {correct}");
		recip(temp, alp, alpinv, temp2);
		Console.WriteLine($"Step P done at correct = {correct}");
		mpmult(temp2, temp3, alp, alpinv, pk1);
		Console.WriteLine($"Step Q done at correct = {correct}");

		/*
		This for loop prints out pi so that the user can see the
		convergence.
		*/

		Console.WriteLine($"This is pi to the desired precision at correct = {correct}");
		Console.WriteLine();

		for (i = (pk1[5002] - 1); i > -1; i--)
		{
			if (pk1[i] >= 10)
			{
				Console.Write(pk1[i]);
			}
			else
			{
				Console.Write($"0{pk1[i]}");
			}
		}

		Console.WriteLine();

		/*
		This for loop sets the i+1 element as i.
		*/

		for (i = 0; i < 5003; i++)
		{
			ak[i] = ak1[i];
			bk[i] = bk1[i];
			pk[i] = pk1[i];
		}

		for (i = 0; i < 5003; i++)
		{
			temp[i] = 0;
			temp2[i] = 0;
			temp3[i] = 0;
			ak1[i] = 0;
			bk1[i] = 0;
			pk1[i] = 0;
		}

		if (correct >= 1)
		{
			correct = (short)(correct * 2);
		}
		else
		{
			correct = 1;
		}
	}

	/*
	This for loop print out pi to the desired precision. 
	*/

	Console.WriteLine("This is pi to the desired precision");
	Console.WriteLine();

	var piBuilder = new StringBuilder();

	for (i = (pk[5002] - 1); i > -1; i--)
	{
		if (pk[i] > 10)
		{
			Console.Write(pk[i]);
			piBuilder.Append(pk[i]);
		}
		else
		{
			Console.Write($"0{pk[i]}");
			piBuilder.Append(CultureInfo.CurrentCulture, $"0{pk[i]}");
		}
	}

	File.WriteAllText("pi.txt", piBuilder.ToString());
}

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
	int n0, n1, n2, n3, n4, n5, n6, k0, k1, k2, k3, k4, k5, k6, accout, accin, accrt;
	int t, temp;
	int[] x0, x1, h0, h1, g0;

	long check, go = int.MaxValue;
	short i;

	x0 = new int[5003];
	x1 = new int[5003];
	h0 = new int[5003];
	h1 = new int[5003];
	g0 = new int[5003];

	for (i = 0; i < 5003; i++)
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
										check = Math.Abs(go / alp[accrt]);

										if (Math.Abs(vect1[accin]) > check)
										{
											temp = mymod(vect1[accin], alp[accrt]);
											x1[accout] = negmod(x1[accout] + temp, M);
										}
										else
										{
											temp = negmod(vect1[accin] * alp[accrt], M);
											x1[accout] = negmod(x1[accout] + temp, M);
										}

										if (Math.Abs(vect2[accin]) > check)
										{
											temp = mymod(vect2[accin], alp[accrt]);
											h1[accout] = negmod(h1[accout] + temp, M);
										}
										else
										{
											temp = negmod(vect2[accin] * alp[accrt], M);
											h1[accout] = negmod(h1[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		h0[i] = 0;
		x0[i] = 0;
	}

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
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k5 = 0; k5 < 2; k5++)
								{
									accout = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * k3 + 4 * k4 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * k3 + 4 * k4 + 2 * k5 + n0;
									accrt = (1250 * k5 * (2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										x0[accout] = negmod(x1[accin] + x0[accout], M);
										h0[accout] = negmod(h1[accin] + h0[accout], M);
									}
									else if (accrt == 2500)
									{
										x0[accout] = negmod(-x1[accin] + x0[accout], M);
										h0[accout] = negmod(-h1[accin] + h0[accout], M);
									}
									else
									{
										check = Math.Abs(go / alp[accrt]);

										if (Math.Abs(x1[accin]) > check)
										{
											temp = mymod(x1[accin], alp[accrt]);
											x0[accout] = negmod(x0[accout] + temp, M);
										}
										else
										{
											temp = negmod(x1[accin] * alp[accrt], M);
											x0[accout] = negmod(x0[accout] + temp, M);
										}

										if (Math.Abs(h1[accin]) > check)
										{
											temp = mymod(h1[accin], alp[accrt]);
											h0[accout] = negmod(h0[accout] + temp, M);
										}
										else
										{
											temp = negmod(h1[accin] * alp[accrt], M);
											h0[accout] = negmod(h0[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		x1[i] = 0;
		h1[i] = 0;
	}

	for (k0 = 0; k0 < 5; k0++)
	{
		for (k1 = 0; k1 < 5; k1++)
		{
			for (k2 = 0; k2 < 5; k2++)
			{
				for (k3 = 0; k3 < 5; k3++)
				{
					for (n2 = 0; n2 < 2; n2++)
					{
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k4 = 0; k4 < 2; k4++)
								{
									accout = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * k3 + 4 * n2 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * k3 + 4 * k4 + 2 * n1 + n0;
									accrt = (625 * k4 * (4 * n2 + 2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										x1[accout] = negmod(x0[accin] + x1[accout], M);
										h1[accout] = negmod(h0[accin] + h1[accout], M);
									}
									else if (accrt == 2500)
									{
										x1[accout] = negmod(-x0[accin] + x1[accout], M);
										h1[accout] = negmod(-h0[accin] + h1[accout], M);
									}
									else
									{
										check = Math.Abs(go / alp[accrt]);

										if (Math.Abs(x0[accin]) > check)
										{
											temp = mymod(x0[accin], alp[accrt]);
											x1[accout] = negmod(x1[accout] + temp, M);
										}
										else
										{
											temp = negmod(x0[accin] * alp[accrt], M);
											x1[accout] = negmod(x1[accout] + temp, M);
										}

										if (Math.Abs(h0[accin]) > check)
										{
											temp = mymod(h0[accin], alp[accrt]);
											h1[accout] = negmod(h1[accout] + temp, M);
										}
										else
										{
											temp = negmod(h0[accin] * alp[accrt], M);
											h1[accout] = negmod(h1[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		x0[i] = 0;
		h0[i] = 0;
	}

	for (k0 = 0; k0 < 5; k0++)
	{
		for (k1 = 0; k1 < 5; k1++)
		{
			for (k2 = 0; k2 < 5; k2++)
			{
				for (n3 = 0; n3 < 5; n3++)
				{
					for (n2 = 0; n2 < 2; n2++)
					{
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k3 = 0; k3 < 5; k3++)
								{
									accout = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * k3 + 4 * n2 + 2 * n1 + n0;
									accrt = (125 * k3 * (8 * n3 + 4 * n2 + 2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										x0[accout] = negmod(x1[accin] + x0[accout], M);
										h0[accout] = negmod(h1[accin] + h0[accout], M);
									}
									else if (accrt == 2500)
									{
										x0[accout] = negmod(-x1[accin] + x0[accout], M);
										h0[accout] = negmod(-h1[accin] + h0[accout], M);
									}
									else
									{
										check = Math.Abs(go / alp[accrt]);

										if (Math.Abs(x1[accin]) > check)
										{
											temp = mymod(x1[accin], alp[accrt]);
											x0[accout] = negmod(x0[accout] + temp, M);
										}
										else
										{
											temp = negmod(x1[accin] * alp[accrt], M);
											x0[accout] = negmod(x0[accout] + temp, M);
										}

										if (Math.Abs(h1[accin]) > check)
										{
											temp = mymod(h1[accin], alp[accrt]);
											h0[accout] = negmod(h0[accout] + temp, M);
										}
										else
										{
											temp = negmod(h1[accin] * alp[accrt], M);
											h0[accout] = negmod(h0[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		x1[i] = 0;
		h1[i] = 0;
	}

	for (k0 = 0; k0 < 5; k0++)
	{
		for (k1 = 0; k1 < 5; k1++)
		{
			for (n4 = 0; n4 < 5; n4++)
			{
				for (n3 = 0; n3 < 5; n3++)
				{
					for (n2 = 0; n2 < 2; n2++)
					{
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k2 = 0; k2 < 5; k2++)
								{
									accout = 1000 * k0 + 200 * k1 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accrt = (25 * k2 * (40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										x1[accout] = negmod(x0[accin] + x1[accout], M);
										h1[accout] = negmod(h0[accin] + h1[accout], M);
									}
									else if (accrt == 2500)
									{
										x1[accout] = negmod(-x0[accin] + x1[accout], M);
										h1[accout] = negmod(-h0[accin] + h1[accout], M);
									}
									else
									{
										check = Math.Abs(go / alp[accrt]);

										if (Math.Abs(x0[accin]) > check)
										{
											temp = mymod(x0[accin], alp[accrt]);
											x1[accout] = negmod(x1[accout] + temp, M);
										}
										else
										{
											temp = negmod(x0[accin] * alp[accrt], M);
											x1[accout] = negmod(x1[accout] + temp, M);
										}

										if (Math.Abs(h0[accin]) > check)
										{
											temp = mymod(h0[accin], alp[accrt]);
											h1[accout] = negmod(h1[accout] + temp, M);
										}
										else
										{
											temp = negmod(h0[accin] * alp[accrt], M);
											h1[accout] = negmod(h1[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		x0[i] = 0;
		h0[i] = 0;
	}

	for (k0 = 0; k0 < 5; k0++)
	{
		for (n5 = 0; n5 < 5; n5++)
		{
			for (n4 = 0; n4 < 5; n4++)
			{
				for (n3 = 0; n3 < 5; n3++)
				{
					for (n2 = 0; n2 < 2; n2++)
					{
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k1 = 0; k1 < 5; k1++)
								{
									accout = 1000 * k0 + 200 * n5 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * k1 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accrt = (5 * k1 * (200 * n5 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										x0[accout] = negmod(x1[accin] + x0[accout], M);
										h0[accout] = negmod(h1[accin] + h0[accout], M);
									}
									else if (accrt == 2500)
									{
										x0[accout] = negmod(-x1[accin] + x0[accout], M);
										h0[accout] = negmod(-h1[accin] + h0[accout], M);
									}
									else
									{
										check = Math.Abs(go / alp[accrt]);

										if (Math.Abs(x1[accin]) > check)
										{
											temp = mymod(x1[accin], alp[accrt]);
											x0[accout] = negmod(x0[accout] + temp, M);
										}
										else
										{
											temp = negmod(x1[accin] * alp[accrt], M);
											x0[accout] = negmod(x0[accout] + temp, M);
										}

										if (Math.Abs(h1[accin]) > check)
										{
											temp = mymod(h1[accin], alp[accrt]);
											h0[accout] = negmod(h0[accout] + temp, M);
										}
										else
										{
											temp = negmod(h1[accin] * alp[accrt], M);
											h0[accout] = negmod(h0[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		x1[i] = 0;
		h1[i] = 0;
	}

	for (n6 = 0; n6 < 5; n6++)
	{
		for (n5 = 0; n5 < 5; n5++)
		{
			for (n4 = 0; n4 < 5; n4++)
			{
				for (n3 = 0; n3 < 5; n3++)
				{
					for (n2 = 0; n2 < 2; n2++)
					{
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k0 = 0; k0 < 5; k0++)
								{
									accout = 1000 * n6 + 200 * n5 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * n5 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accrt = (k0 * (1000 * n6 + 200 * n5 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										x1[accout] = negmod(x0[accin] + x1[accout], M);
										h1[accout] = negmod(h0[accin] + h1[accout], M);
									}
									else if (accrt == 2500)
									{
										x1[accout] = negmod(-x0[accin] + x1[accout], M);
										h1[accout] = negmod(-h0[accin] + h1[accout], M);
									}
									else
									{
										check = Math.Abs(go / alp[accrt]);

										if (Math.Abs(x0[accin]) > check)
										{
											temp = mymod(x0[accin], alp[accrt]);
											x1[accout] = negmod(x1[accout] + temp, M);
										}
										else
										{
											temp = negmod(x0[accin] * alp[accrt], M);
											x1[accout] = negmod(x1[accout] + temp, M);
										}

										if (Math.Abs(h0[accin]) > check)
										{
											temp = mymod(h0[accin], alp[accrt]);
											h1[accout] = negmod(h1[accout] + temp, M);
										}
										else
										{
											temp = negmod(h0[accin] * alp[accrt], M);
											h1[accout] = negmod(h1[accout] + temp, M);
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

	/*
	This for loop multiplies the two transformed vectors together element by
	element
	*/

	for (i = 0; i < 5000; i++)
	{
		if (h1[i] != 0)
		{
			check = Math.Abs(go / h1[i]);

			if (Math.Abs(x1[i]) > check)
			{
				g0[i] = mymod(x1[i], h1[i]);
			}
			else
			{
				g0[i] = negmod(x1[i] * h1[i], M);
			}
		}
		else
		{
			g0[i] = 0;
		}
	}

	/*
	So far, x1 and h1 holds the transform of the input vector x0; also,
	x0[5000 -> 5002] and h0[5000 -> 5002] have not been touched and will not
	for the remainder of the calculation, which performs the inverse transform of
	g0[].
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
										g1[accout] = negmod(g0[accin] + g1[accout], M);
									}
									else if (accrt == 2500)
									{
										g1[accout] = negmod(-g0[accin] + g1[accout], M);
									}
									else
									{
										check = Math.Abs(go / alpinv[accrt]);

										if (Math.Abs(g0[accin]) > check)
										{
											temp = mymod(g0[accin], alpinv[accrt]);
											g1[accout] = negmod(g1[accout] + temp, M);
										}
										else
										{
											temp = negmod(g0[accin] * alpinv[accrt], M);
											g1[accout] = negmod(g1[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		g0[i] = 0;
	}

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
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k5 = 0; k5 < 2; k5++)
								{
									accout = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * k3 + 4 * k4 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * k3 + 4 * k4 + 2 * k5 + n0;
									accrt = (1250 * k5 * (2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										g0[accout] = negmod(g1[accin] + g0[accout], M);
									}
									else if (accrt == 2500)
									{
										g0[accout] = negmod(-g1[accin] + g0[accout], M);
									}
									else
									{
										check = Math.Abs(go / alpinv[accrt]);

										if (Math.Abs(g1[accin]) > check)
										{
											temp = mymod(g1[accin], alpinv[accrt]);
											g0[accout] = negmod(g0[accout] + temp, M);
										}
										else
										{
											temp = negmod(g1[accin] * alpinv[accrt], M);
											g0[accout] = negmod(g0[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		g1[i] = 0;
	}

	for (k0 = 0; k0 < 5; k0++)
	{
		for (k1 = 0; k1 < 5; k1++)
		{
			for (k2 = 0; k2 < 5; k2++)
			{
				for (k3 = 0; k3 < 5; k3++)
				{
					for (n2 = 0; n2 < 2; n2++)
					{
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k4 = 0; k4 < 2; k4++)
								{
									accout = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * k3 + 4 * n2 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * k3 + 4 * k4 + 2 * n1 + n0;
									accrt = (625 * k4 * (4 * n2 + 2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										g1[accout] = negmod(g0[accin] + g1[accout], M);
									}
									else if (accrt == 2500)
									{
										g1[accout] = negmod(-g0[accin] + g1[accout], M);
									}
									else
									{
										check = Math.Abs(go / alpinv[accrt]);

										if (Math.Abs(g0[accin]) > check)
										{
											temp = mymod(g0[accin], alpinv[accrt]);
											g1[accout] = negmod(g1[accout] + temp, M);
										}
										else
										{
											temp = negmod(g0[accin] * alpinv[accrt], M);
											g1[accout] = negmod(g1[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		g0[i] = 0;
	}

	for (k0 = 0; k0 < 5; k0++)
	{
		for (k1 = 0; k1 < 5; k1++)
		{
			for (k2 = 0; k2 < 5; k2++)
			{
				for (n3 = 0; n3 < 5; n3++)
				{
					for (n2 = 0; n2 < 2; n2++)
					{
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k3 = 0; k3 < 5; k3++)
								{
									accout = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * k3 + 4 * n2 + 2 * n1 + n0;
									accrt = (125 * k3 * (8 * n3 + 4 * n2 + 2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										g0[accout] = negmod(g1[accin] + g0[accout], M);
									}
									else if (accrt == 2500)
									{
										g0[accout] = negmod(-g1[accin] + g0[accout], M);
									}
									else
									{
										check = Math.Abs(go / alpinv[accrt]);

										if (Math.Abs(g1[accin]) > check)
										{
											temp = mymod(g1[accin], alpinv[accrt]);
											g0[accout] = negmod(g0[accout] + temp, M);
										}
										else
										{
											temp = negmod(g1[accin] * alpinv[accrt], M);
											g0[accout] = negmod(g0[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		g1[i] = 0;
	}

	for (k0 = 0; k0 < 5; k0++)
	{
		for (k1 = 0; k1 < 5; k1++)
		{
			for (n4 = 0; n4 < 5; n4++)
			{
				for (n3 = 0; n3 < 5; n3++)
				{
					for (n2 = 0; n2 < 2; n2++)
					{
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k2 = 0; k2 < 5; k2++)
								{
									accout = 1000 * k0 + 200 * k1 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * k1 + 40 * k2 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accrt = (25 * k2 * (40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										g1[accout] = negmod(g0[accin] + g1[accout], M);
									}
									else if (accrt == 2500)
									{
										g1[accout] = negmod(-g0[accin] + g1[accout], M);
									}
									else
									{
										check = Math.Abs(go / alpinv[accrt]);

										if (Math.Abs(g0[accin]) > check)
										{
											temp = mymod(g0[accin], alpinv[accrt]);
											g1[accout] = negmod(g1[accout] + temp, M);
										}
										else
										{
											temp = negmod(g0[accin] * alpinv[accrt], M);
											g1[accout] = negmod(g1[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		g0[i] = 0;
	}

	for (k0 = 0; k0 < 5; k0++)
	{
		for (n5 = 0; n5 < 5; n5++)
		{
			for (n4 = 0; n4 < 5; n4++)
			{
				for (n3 = 0; n3 < 5; n3++)
				{
					for (n2 = 0; n2 < 2; n2++)
					{
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k1 = 0; k1 < 5; k1++)
								{
									accout = 1000 * k0 + 200 * n5 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * k1 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accrt = (5 * k1 * (200 * n5 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										g0[accout] = negmod(g1[accin] + g0[accout], M);
									}
									else if (accrt == 2500)
									{
										g0[accout] = negmod(-g1[accin] + g0[accout], M);
									}
									else
									{
										check = Math.Abs(go / alpinv[accrt]);

										if (Math.Abs(g1[accin]) > check)
										{
											temp = mymod(g1[accin], alpinv[accrt]);
											g0[accout] = negmod(g0[accout] + temp, M);
										}
										else
										{
											temp = negmod(g1[accin] * alpinv[accrt], M);
											g0[accout] = negmod(g0[accout] + temp, M);
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

	for (i = 0; i < 5000; i++)
	{
		g1[i] = 0;
	}

	for (n6 = 0; n6 < 5; n6++)
	{
		for (n5 = 0; n5 < 5; n5++)
		{
			for (n4 = 0; n4 < 5; n4++)
			{
				for (n3 = 0; n3 < 5; n3++)
				{
					for (n2 = 0; n2 < 2; n2++)
					{
						for (n1 = 0; n1 < 2; n1++)
						{
							for (n0 = 0; n0 < 2; n0++)
							{
								for (k0 = 0; k0 < 5; k0++)
								{
									accout = 1000 * n6 + 200 * n5 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accin = 1000 * k0 + 200 * n5 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0;
									accrt = (k0 * (1000 * n6 + 200 * n5 + 40 * n4 + 8 * n3 + 4 * n2 + 2 * n1 + n0)) % N;

									if (accrt == 0)
									{
										g1[accout] = negmod(g0[accin] + g1[accout], M);
									}
									else if (accrt == 2500)
									{
										g1[accout] = negmod(-g0[accin] + g1[accout], M);
									}
									else
									{
										check = Math.Abs(go / alpinv[accrt]);

										if (Math.Abs(g0[accin]) > check)
										{
											temp = mymod(g0[accin], alpinv[accrt]);
											g1[accout] = negmod(g1[accout] + temp, M);
										}
										else
										{
											temp = negmod(g0[accin] * alpinv[accrt], M);
											g1[accout] = negmod(g1[accout] + temp, M);
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

	/*
	N-1 must be multiplied to each element in g1[]
	*/

	for (i = 0; i < 5000; i++)
	{
		check = Math.Abs(go / -10002);

		if (Math.Abs(g1[i]) > check)
		{
			g1[i] = mymod(g1[i], -10002);
		}
		else
		{
			g1[i] = negmod(g1[i] * -10002, M);
		}
	}

	/*
	This for loop makes all of the elements in g1[] positive. 
	*/

	for (i = 0; i < 5000; i++)
	{
		if (g1[i] < 0)
		{
			g1[i] = g1[i] + M;
		}
	}

	/*
	This for loop released the carries on g1[] and also determines the
	exponent and length of g1[].
	*/

	for (i = 0; i < (vect1[5002] + vect2[5002]); i++)
	{
		t = g1[i] / 100;
		temp = t * 100;
		g1[i] = g1[i] - temp;
		g1[i + 1] = g1[i + 1] + t;

		if (i == (vect1[5002] + vect2[5002] - 2) || (vect1[5002] + vect2[5002]) < 2)
		{
			if (t == 0)
			{
				g1[5001] = vect1[5001] + vect2[5001];

				if ((vect1[5002] + vect2[5002]) < 2)
				{
					g1[5002] = 1;
				}
				else
				{
					g1[5002] = vect1[5002] + vect2[5002] - 1;
				}
			}
			else
			{
				g1[5001] = vect1[5001] + vect2[5001] + 1;

				if ((vect1[5002] + vect2[5002]) < 2)
				{
					g1[5002] = 2;
				}
				else
				{
					g1[5002] = vect1[5002] + vect2[5002];
				}
			}
		}
	}

	/*
	This for loop sets all elements that do not contain part of the answer
	to 0.
	*/

	for (i = (short)g1[5002]; i < 5000; i++)
	{
		g1[i] = 0;
	}

	for (i = 0; i < 5003; i++)
	{
		x0[i] = 0;
	}

	/*
	This for loop truncates the final vector to the limit value defined by
	the user if the answer in g1[] is larger than limit.
	*/

	if (g1[5002] > Shared.limit)
	{
		temp = g1[5002] - Shared.limit;

		for (i = (short)((short)(g1[5002]) - 1); i > temp - 1; i--)
		{
			x0[i - temp] = g1[i];
		}

		for (i = 0; i < 5000; i++)
		{
			g1[i] = x0[i];
		}

		g1[5002] = Shared.limit;
	}

	/*
	This sets the sign on g1[] to be either positive or negative;
	*/

	if (vect1[5000] == vect2[5000])
	{
		g1[5000] = 1;
	}
	else
	{
		g1[5000] = 0;
	}
}

static void mpadd(int[] vect1, int[] vect2, int[] g)
{
	short i;
	int t, temp, shift, tell, length, maxexp;
	var x = new int[5003];
	var h = new int[5003];

	/*
	Clear out the vectors;
	*/

	for (i = 0; i < 5003; i++)
	{
		g[i] = 0;
		x[i] = 0;
		h[i] = 0;
	}

	maxexp = Math.Max(vect1[5001], vect2[5001]);

	/*
	This if-else if-else structure determines how the input vectors should
	be shifted (if at all) to align the numbers such that the decimal point
	is located in the same place in the vectors.
	*/

	if ((vect1[5002] - vect1[5001]) < (vect2[5002] - vect2[5001]))
	{
		shift = vect2[5002] - vect1[5002] + vect1[5001] - vect2[5001];

		if ((vect1[5002] + shift) >= vect2[5002])
		{
			length = vect1[5002] + shift;
		}
		else
		{
			length = vect2[5002];
		}

		tell = 1;
	}
	else if ((vect1[5002] - vect1[5001]) > (vect2[5002] - vect2[5001]))
	{
		shift = vect1[5002] - vect2[5002] + vect2[5001] - vect1[5001];

		if ((vect2[5002] + shift) >= vect1[5002])
		{
			length = vect2[5002] + shift;
		}
		else
		{
			length = vect1[5002];
		}

		tell = 0;
	}
	else
	{
		shift = 0;
		tell = 2;
		length = Math.Max(vect1[5002], vect2[5002]);
	}

	/*
	These if statements shift the input vectors (if necessary). 
	*/

	if (tell == 1)
	{
		for (i = 0; i < vect1[5002]; i++)
		{
			x[i + shift] = vect1[i];
		}
	}
	else
	{
		for (i = 0; i < vect1[5002]; i++)
		{
			x[i] = vect1[i];
		}
	}

	if (tell == 0)
	{
		for (i = 0; i < vect2[5002]; i++)
		{
			h[i + shift] = vect2[i];
		}
	}
	else
	{
		for (i = 0; i < vect2[5002]; i++)
		{
			h[i] = vect2[i];
		}
	}

	/*
	This for loop adds the x and h vectors together.
	*/

	for (i = 0; i < length; i++)
	{
		g[i] = x[i] + h[i];
	}

	/*
	This for loop released the carries on g[].
	*/

	for (i = 0; i < length; i++)
	{
		t = g[i] / 100;
		temp = t * 100;
		g[i] = g[i] - temp;
		g[i + 1] = g[i + 1] + t;

		if (i == (length - 2) || length < 2)
		{
			if (t == 0)
			{
				g[5001] = maxexp;
				g[5002] = length;
			}
			else
			{
				g[5001] = maxexp + 1;
				g[5002] = length + 1;
			}
		}
	}

	/*
	This for loop sets all elements that do not contain a part of the answer
	to 0.
	*/

	for (i = (short)g[5002]; i < 5000; i++)
	{
		g[i] = 0;
	}

	for (i = 0; i < 5003; i++)
	{
		x[i] = 0;
	}

	/*
	This for loop truncates the final vector to the limit value defined by
	the user if the answer in g1[] is larger than limit.
	*/

	if (g[5002] > Shared.limit)
	{
		temp = g[5002] - Shared.limit;

		for (i = (short)(g[5002] - 1); i > temp - 1; i--)
		{
			x[i - temp] = g[i];

			for (i = 0; i < 5000; i++)
			{
				g[i] = x[i];
			}

			g[5002] = Shared.limit;
		}
	}

	g[5000] = 1;
}

/*
The mpsub function subtracts vect2[] from vect1[]
*/

static void mpsub(int[] vect1, int[] vect2, int[] g)
{
	short i, j, set;
	int stemp, temp, shift, tell, length, maxexp;

	var x = new int[5003];
	var h = new int[5003];

	set = 0;

	/*
	Clear out the vectors 
	*/

	for (i = 0; i < 5003; i++)
	{
		x[i] = 0;
		h[i] = 0;
		g[i] = 0;
	}

	maxexp = Math.Max(vect1[5001], vect2[5001]);

	/*
	This if-else if-else structure determines how th einput vectors should
	be shifted (if at all) to align the numbers such that the decimal point
	is located in the same place in the vectors.
	*/

	if ((vect1[5002] - vect1[5001]) < (vect2[5002 - vect2[5001]]))
	{
		shift = vect2[5002] - vect1[5002] + vect1[5001] - vect2[5001];

		if ((vect1[5002] + shift) >= vect2[5002])
		{
			length = vect1[5002] + shift;
		}
		else
		{
			length = vect2[5002];
		}

		tell = 1;
	}
	else if ((vect1[5002] - vect1[5001]) > (vect2[5002] - vect2[5001]))
	{
		shift = vect1[5002] - vect2[5002] + vect2[5001] - vect1[5001];

		if ((vect2[5002] + shift) >= vect1[5002])
		{
			length = vect2[5002] + shift;
		}
		else
		{
			length = vect2[5002];
		}

		tell = 0;
	}
	else
	{
		shift = 0;
		tell = 2;
		length = Math.Max(vect1[5002], vect2[5002]);
	}

	/*
	These if statements shift the input vectors (if necessary). 
	*/

	if (tell == 1)
	{
		for (i = 0; i < vect1[5002]; i++)
		{
			x[i + shift] = vect1[i];
		}
	}
	else
	{
		for (i = 0; i < vect1[5002]; i++)
		{
			x[i] = vect1[i];
		}
	}

	if (tell == 0)
	{
		for (i = 0; i < vect2[5002]; i++)
		{
			h[i + shift] = vect2[i];
		}
	}
	else
	{
		for (i = 0; i < vect2[5002]; i++)
		{
			h[i] = vect2[i];
		}
	}

	/*
	This while loop determines if vect2[] was bigger than vect1[]. If so
	then the result will be negative; otherwise, the answer will be positive.
	If the result will be negative, the x and h vectors are switched.
	*/

	i = (short)(length - 1);

	while (i > -1)
	{
		temp = x[i] - h[1];

		if (temp > 0)
		{
			g[5000] = 1;
			i = -2;
		}
		else if (temp < 0)
		{
			for (j = (short)length; j > -1; j--)
			{
				stemp = x[j];
				x[j] = h[j];
				h[j] = stemp;
				g[5000] = 0;
				i = -2;
			}
		}
		else
		{
			i--;
		}
	}

	/*
	This for loop subtracts h[] from x[] and stores the result in g[].
	*/

	for (i = 0; i < length; i++)
	{
		if (x[i] >= h[i])
		{
			g[i] = x[i] - h[i];
		}
		else
		{
			x[i + 1] = x[i + 1] - 1;
			x[i] = x[i] + 100;
			g[i] = x[i] - h[i];
		}
	}

	/*
	This while loop determines the exponent and length of g[].
	*/

	i = (short)(length - 1);

	while (i > -1)
	{
		if (g[i] != 0)
		{
			i = 2;
		}
		else
		{
			if (length == 1)
			{
				maxexp = 0;
				length = 1;
				i--;
				set = 1;
			}
			else
			{
				maxexp = maxexp - 1;
				length = length - 1;
				i--;
			}
		}
	}

	if (set == 1)
	{
		g[5000] = 1;
	}

	g[5001] = maxexp;
	g[5002] = length;

	/*
	This for loop sets all elements that do not contain part of the answer
	to 0.
	*/

	for (i = (short)g[5002]; i < 5000; i++)
	{
		g[i] = 0;
	}

	for (i = 0; i < 5003; i++)
	{
		x[i] = 0;
	}

	/*
	This for loop truncates the final vector to the limit value defined by
	the user if the answer in g1[] is larger than limit.
	*/

	if (g[5002] > Shared.limit)
	{
		temp = g[5002] - Shared.limit;

		for (i = (short)(g[5002] - 1); i > (temp - 1); i--)
		{
			x[i - temp] = g[i];
		}

		for (i = 0; i < 5000; i++)
		{
			g[i] = x[i];
		}

		g[5002] = Shared.limit;
	}
}

/*
This function calculates an initial guess of the square root of
vect[].
*/

static void sqrtguess(int[] vect, int[] g)
{
	short i;
	int k, yguess;
	long temp, one;

	/*
	Clear out the vector.
	*/
	for (i = 0; i < 5003; i++)
	{
		g[i] = 0;
	}

	/*
	This if structure calculates the initial guess of vect1[].
	*/

	if (vect[5001] % 2 == 0 && vect[5002] >= 3)
	{
		k = vect[5001] / 2;
		temp = (10000 * vect[vect[5002] - 1] + 100 * vect[vect[5002] - 2] + vect[vect[5002] - 3]);
		yguess = (int)Math.Sqrt(temp);
		g[0] = yguess - (yguess / 100) * 100;
		yguess = yguess / 100;
		g[1] = yguess - (yguess / 100) * 100;
		g[2] = yguess / 100;
		g[5000] = 1;
		g[5001] = k;

		if (g[2] != 0)
		{
			g[5002] = 3;
		}
		else
		{
			g[5002] = 2;
		}
	}
	else if (vect[5001] % 2 != 0 && vect[5002] >= 2)
	{
		k = vect[5001] - 1;
		k = k / 2;
		temp = (100 * vect[vect[5002] - 1] + vect[vect[5002] - 2]);
		yguess = (int)Math.Sqrt(temp);
		g[0] = yguess - (yguess / 100) * 100;
		g[1] = yguess / 100;
		g[5000] = 1;
		g[5001] = k;

		if (g[1] != 0)
		{
			g[5002] = 2;
		}
		else
		{
			g[5002] = 1;
		}
	}
	else if (vect[5002] == 2)
	{
		if (vect[5001] % 2 == 0)
		{
			k = vect[5001] / 2;
		}
		else
		{
			k = vect[5001] - 1;
			k = k / 2;
		}

		temp = (100 * vect[1] + vect[0]);
		yguess = (int)Math.Sqrt(temp);
		g[0] = yguess - (yguess / 100) * 100;
		g[1] = yguess / 100;
		g[5000] = 1;
		g[5001] = k;

		if (g[1] != 0)
		{
			g[5002] = 2;
		}
		else
		{
			g[5002] = 1;
		}
	}
	else
	{
		if (vect[5001] % 2 == 0)
		{
			k = vect[5001] / 2;
		}
		else
		{
			k = vect[5001] - 1;
			k = k / 2;
		}

		one = (long)Math.Sqrt(vect[0]);
		g[1] = (int)one;
		g[0] = (int)((one - g[1]) * 100);
		g[5000] = 1;
		g[5001] = k;

		if (g[1] != 0)
		{
			g[5002] = 2;
		}
		else
		{
			g[5002] = 1;
		}
	}
}

/*
The ptfive function multiplies in[] by 0.5. Since the recipsqrt
function as well as the main program needs to divide a vector by two,
this function was written to perform the multiplication rather than 
calling mpmult (this is quicker).
*/

static void ptfive(int[] @in, int[] @out)
{
	short i;
	int t, temp;

	/*
	Clear out the vector.
	*/

	for (i = 0; i < 5003; i++)
	{
		@out[i] = 0;
	}

	/*
	This for loop performs the multiplication.
	*/

	for (i = 0; i < @in[5002]; i++)
	{
		@out[i] = @in[i] * 50;
	}

	/*
	This for loop releases the carries on @out[].
	*/

	for (i = 0; i < @in[5002]; i++)
	{
		t = @out[i] / 100;
		temp = t * 100;
		@out[i] = @out[i] - temp;
		@out[i + 1] = @out[i + 1] + t;

		if (i == (@in[5002] - 1) || @in[5002] < 2)
		{
			if (t == 0)
			{
				@out[5001] = @in[5001] = 1;

				if (@in[5002] < 2)
				{
					@out[5002] = 1;
				}
				else
				{
					@out[5002] = @in[5002];
				}
			}
			else
			{
				@out[5001] = @in[5001];

				if (@in[5002] < 2)
				{
					@out[5002] = 2;
				}
				else
				{
					@out[5002] = @in[5002] + 1;
				}
			}
		}
	}

	if (@in[5000] == 1)
	{
		@out[5000] = 1;
	}
	else
	{
		@out[5000] = 0;
	}
}

/*
The invert function performs 1 / vect[] for the recipsqrt function.
*/
static void invert(int[] vect, int[] g)
{
	short i;
	double show;

	/*
	Clear out the vector.
	*/

	for (i = 0; i < 5003; i++)
	{
		g[i] = 0;
	}

	show = 100 * vect[1] + vect[0];

	if (show == 100 && vect[5001] == 0)
	{
		g[0] = 1;
		g[5002] = 1;
		g[5000] = 1;
		g[5001] = 0;
	}
	else
	{
		show = 1 / show;

		g[1] = (int)(100 * show);

		if (g[1] == 0)
		{
			show = show * 100;
			g[1] = (int)(100 * show);
		}

		show = show * 100 - 100 * show;
		g[0] = (int)(100 * show);
		g[5002] = 2;
		g[5000] = 1;
		g[5001] = vect[5001] * -1 - 1;
	}
}

/*
The recipguess function performs 1/vect[] for the recip function.
*/

static void recipguess(int[] vect, int[] g)
{
	short i;
	double show;

	/*
	Clear out the vector.
	*/

	for (i = 0; i < 5003; i++)
	{
		g[i] = 0;
	}

	if (vect[5002] > 1)
	{
		show = 100.0 * vect[vect[5002] - 1] + vect[vect[5002] - 2];
	}
	else
	{
		show = vect[0];
	}

	if (show == 100.0 && vect[5001] == 0)
	{
		g[0] = 1;
		g[5002] = 1;
		g[5000] = 1;
		g[5001] = 0;
	}
	else
	{
		show = 1.0 / show;
		g[1] = (int)(100.0 * show);

		if (g[1] == 0)
		{
			show = show * 100.0;
			g[1] = (int)(100.0 * show);
		}

		show = show * 100.0 - (100.0 * show);
		g[0] = (int)(100.0 * show);
		g[5002] = 2;
		g[5000] = 1;
		g[5001] = vect[5001] * (-1) - 1;
	}
}

/*
The recip function calculates 1/y[] using Newton's Method.
*/

static void recip(int[] y, int[] alp, int[] alpinv, int[] first)
{
	var two = new int[5003];
	var second = new int[5003];
	var third = new int[5003];
	var fourth = new int[5003];

	short correct = 0, i;

	/*
	Clear out the vectors.
	*/
	for (i = 0; i < 5003; i++)
	{
		two[i] = 0;
		first[i] = 0;
		second[i] = 0;
		third[i] = 0;
		fourth[i] = 0;
	}

	/*
	Define vector two
	*/

	two[0] = 2;
	two[5000] = 1;
	two[5001] = 0;
	two[5002] = 1;

	recipguess(y, first);

	while (correct < Shared.limit)
	{
		mpmult(first, y, alp, alpinv, second);

		if (second[5000] == 1)
		{
			mpsub(two, second, third);
		}
		else
		{
			second[5000] = 1;
			mpadd(two, second, third);
		}

		mpmult(third, first, alp, alpinv, fourth);

		for (i = 0; i < 5003; i++)
		{
			first[i] = fourth[i];
		}

		for (i = 0; i < 5003; i++)
		{
			second[i] = 0;
			third[i] = 0;
			fourth[i] = 0;
		}

		if (correct >= 1)
		{
			correct = (short)(correct * 2);
		}
		else
		{
			correct = 1;
		}
	}
}

/*
The recipsqrt function calculates 1/sqrt(y[]) using Newton's Method.
*/

static void recipsqrt(int[] y, int[] alp, int[] alpinv, int[] first)
{
	int i;
	var three = new int[5003];
	var second = new int[5003];
	var third = new int[5003];
	var fourth = new int[5003];
	var fifth = new int[5003];
	var sixth = new int[5003];
	var @try = new int[5003];

	short correct = 0;

	/*
	Clear out the vectors.
	*/

	for (i = 0; i < 5003; i++)
	{
		three[i] = 0;
		first[i] = 0;
		second[i] = 0;
		third[i] = 0;
		fourth[i] = 0;
	}

	/*
	Define vector three 
	*/
	three[0] = 3;
	three[5000] = 1;
	three[5001] = 0;
	three[5002] = 1;

	sqrtguess(y, @try);
	invert(@try, first);

	while (correct < Shared.limit)
	{
		mpmult(first, first, alp, alpinv, second);
		mpmult(second, y, alp, alpinv, third);

		if (second[5000] == 1)
		{
			mpsub(three, third, fourth);
		}
		else
		{
			second[5000] = 1;
			mpadd(three, third, fourth);
		}

		mpmult(first, fourth, alp, alpinv, fifth);
		ptfive(fifth, sixth);

		for (i = 0; i < 5003; i++)
		{
			first[i] = sixth[i];
		}

		for (i = 0; i < 5003; i++)
		{
			second[i] = 0;
			third[i] = 0;
			fourth[i] = 0;
			fifth[i] = 0;
			sixth[i] = 0;
		}

		if (correct >= 1)
		{
			correct = (short)(correct * 2);
		}
		else
		{
			correct = 1;
		}
	}
}

#pragma warning disable CA1050 // Declare types in namespaces
#pragma warning disable CA1716 // Identifiers should not match keywords
public static class Shared
#pragma warning restore CA1716 // Identifiers should not match keywords
#pragma warning restore CA1050 // Declare types in namespaces
{
#pragma warning disable CA1823 // Avoid unused private fields
#pragma warning disable CA2211 // Non-constant fields should not be visible
	public static short limit;
#pragma warning restore CA2211 // Non-constant fields should not be visible
#pragma warning restore CA1823 // Avoid unused private fields
}