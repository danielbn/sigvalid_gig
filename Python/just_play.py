def basic_function(input_a, input_b, input_c, input_d, input_e):

    input_a=input_a + 1
    input_b = input_b + 2
    input_c = input_c + 3
    input_d = input_d + 4
    input_e = input_e + 5

    print 'file processed successfully - GHOST TEST ONLY !!!'
    return (input_a, input_b, input_c, input_d, input_e)
###############

input_a=100
input_b=200
input_c=300
input_d=400
input_e=500

(input_a, input_b, input_c, input_d, input_e)=basic_function(input_a, input_b, input_c, input_d, input_e)

print 'hello world'
print input_c


