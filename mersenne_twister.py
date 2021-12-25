class Mersenne:
    '''
    Made by Levi Verhoef, V2A, 2021
    This class is a translation of the pseudocode at https://en.wikipedia.org/wiki/Mersenne_Twister#Pseudocode
    It uses a mersenne twister to generate random values'''

    def __init__(self, start_seed=5489):
        '''Initialize the coefficients, array, masks and seed'''
        self.w = 32 #wordsize in bits
        self.r = 31 #seperation point of one word
        self.n = 624 #degree of recurrence
        self.m = 397 #middle word, offset
        self.u = 11 #bit shift
        self.s = 7 #tempering bit shift
        self.t = 15 #tempering bit shift
        self.l = 18 #tempering bit shift
        self.a = 0x9908B0DF #coefficient of the rational normal form twist matrix
        self.b = 0x9D2C5680 #tempering bit mask
        self.c = 0xEFC60000 #tempering bit mask
        self.d = 0xFFFFFFFF #additional tempering bit mask
        self.f = 0x6c078965

        self.mersenne_array = [0 for x in range(self.n)]
        self.index = 0

        self.lower_mask = 0x7FFFFFFF
        self.upper_mask = 0x80000000
        self.seed = start_seed
        self.initialize_seed(start_seed)



    def initialize_seed(self, number):
        '''This function initializes the seed'''
        self.mersenne_array[0] = number
        self.index = self.n
        for i in range(1, self.n):
            num = self.f * (self.mersenne_array[i-1] ^ (self.mersenne_array[i-1] >> (self.w-2))) + i
            self.mersenne_array[i] = num & self.d #limit the result to a 32 bit value



    def twist(self):
        '''This function generates the next n values in the series and applies masks to them'''
        for j in range(0, self.n):
            twistnum = (self.mersenne_array[j] & self.upper_mask) + (self.mersenne_array[(j+1) % self.n] & self.lower_mask)
            twistednum = twistnum >> 1
            if (twistnum%2) != 0:
                twistednum = twistednum ^ self.a
            self.mersenne_array[j] = self.mersenne_array[(j + self.m) % self.n] ^ twistednum
        self.index = 0


    def temper_numbers(self):
        '''temper the numbers generated in twist'''
        if self.index >= self.n:
            self.twist()

        y = self.mersenne_array[self.index]
        y = y ^ ((y >> self.u) & self.d)
        y = y ^ ((y << self.s) & self.b)
        y = y ^ ((y << self.t) & self.c)
        y = y ^ (y >> self.l)

        self.index += 1
        return (y & self.d)

    def random(self):
        '''returns a random float value between 0 and 1'''
        return self.temper_numbers() / 2**self.w

    def randint(self, a, b):
        '''returns a random integer between a and b'''
        return int(self.random()/(1/(b-a)) + a)