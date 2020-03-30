"""
Created on March 29, 2020

@author: Yi Wang
"""

def generate_start_end(month):
    """
    month is 'YYYYMM'

    get 'YYYY-MM-01_YYYY-MM-end'
    """

    days_dict = {
            '06' : '30',
            '07' : '31',
            '08' : '31'
            }

    YYYY = month[0:4]
    MM   = month[4:6]

    result = YYYY + '-' + MM + '-01_' + YYYY + '-' + MM + '-' + days_dict[MM]

    return result
