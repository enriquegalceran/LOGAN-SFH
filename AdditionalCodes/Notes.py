# THIS FILE IS NOT TO BE EXECUTED. INSTEAD, IT HAS MULTIPLE SNIPPETS THAT CAN BE USEFUL FOR THE FUTURE


# Argparse that has multiple inputs with variable input parameters.
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', action='append', nargs='+')
args = parser.parse_args()

# In [32]: run test.py -i input1_url input1_name input1_other_var -i input2_url i
# ...: nput2_name input2_other_var -i input3_url input3_name
#
# In [33]: args.i
# Out[33]:
# [['input1_url', 'input1_name', 'input1_other_var'],
#  ['input2_url', 'input2_name', 'input2_other_var'],
#  ['input3_url', 'input3_name']]