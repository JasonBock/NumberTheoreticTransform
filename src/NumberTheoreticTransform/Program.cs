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
	short i, j;
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

// POST-IT, pg. 60, starting with mpsub()...

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















