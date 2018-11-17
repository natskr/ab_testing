
## AB Testing Analysis

* Frequentist AB Testing
* Version df_0: Control
* Version df_1: Contender
* Calculate t-test statistic for continues variables. 
* Calculate z-score for binomial variable.
* Calculate length of test for binomial variables (based on effect_size, type I/II error, test type)
* Calculate P-value: The probability under the null hypothesis of obtaining a result equal to or more extreme than what was observed
* Identify if result is significant


```python
# Read in packages required
import psycopg2
import time
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats import proportion, power

# define significance function
def significant(p, title):
    """Are results significant?"""
    if p <= 0.1: 
          print(title, ': !! Result is Significant !!\n')
    else:
          print(title, ': !! Result is NOT Significant !!\n')

```

### Binomial Variables
* z-score: Deviation from mean in units of standard deviations.
* Like signal-to-noise: (X-m) / SE.
* Sample mean is proportion of success.


```python
# define z-score function for binomial variables
def binomial_var_se(df):
    """Standard Error for binomial variables."""
    binomial_var = float(sum(df.binomial_var_subset) / sum(df.binomial_var_total))
    se = float(np.sqrt(binomial_var*(1-binomial_var) / sum(df.binomial_var_total))) 
    return binomial_var*100, se*100

def z_score_p_value(df, title):
    """z-score for binomial variables."""
    df_0 = df[df.variant_id==0]
    df_1 = df[df.variant_id==1]
    
    binomial_var_0, se_0 = binomial_var_se(df_0)
    binomial_var_1, se_1 = binomial_var_se(df_1)

    z_score = (binomial_var_0 - binomial_var_1) / (np.sqrt((se_0)**2 + (se_1)**2))
    p_values = stats.norm.sf(abs(z_score))*2
    delta = ((binomial_var_1 - binomial_var_0) / binomial_var_0)*100
    output = [title, round(binomial_var_0,5), 
              round(se_0,5), 
              round(binomial_var_1,5), 
              round(se_1,5), 
              df_0.shape[0], 
              df_1.shape[0], 
              round(z_score,5),
              round(p_values,5), 
              round((1-p_values)*100,2), 
              round(delta, 5)]
    print(significant(p_values, title),
          '\ncr_0:', round(binomial_var_0,5), 
          '\nse_0:', round(se_0,5), 
          '\ncr_1:', round(binomial_var_1,5), 
          '\nse_1:', round(se_1,5),
          '\nSize (Control):', df_0.shape[0],
          '\nSize (Contender):', df_1.shape[0],
          '\nz_score:', round(z_score,5),
          '\np_values:', round(p_values,5),
          '\nsignificance:', round((1-p_values)*100,2),
          '\ndelta:', round(delta, 5),'%\n',
          '--------------------------\n')
    return pd.DataFrame([output])


# Calculate length of AB test to be significant
es = proportion.proportion_effectsize(binomial_var_1, 
                                      binomial_var_0,
                                     method='normal')

power.tt_ind_solve_power(effect_size=es,
                        nobs1 = None,
                        alpha=alpha,
                        power=power,
                        ratio=1,
                        alternative='two-sided')

# Calculate P-value
query = """SELECT COUNT(binomial_var_total), 
                  SUM(binomial_var_subset), 
                  variant_id
            FROM table
            GROUP BY variant_id;"""

conn = psycopg2.connect(
        host = host,
        password = password,
        port = port,
        user = user,
        database = database,
    )

cur = conn.cursor()
cur.execute(query)
df = pd.DataFrame(data = cur.fetchall(), 
                  columns = [column[0] for column in cur.description])
conn.close()
print('Number of row:,', df.shape[0])

z_score_p_value(df, 'All')
```

### Continues Variables
* Independent samples (2-tailed) Student's t-test.
* Compares the sample means of 2 groups.


```python
# define t-test for continues variables
def two_sided_t_test(df, title, continues_var, continues_var_name):
    """Two sided t-test for continues variables."""
    df_0 = df[(df.variant_id==0)]
    df_1 = df[(df.variant_id==1)]
    t, p = stats.ttest_ind(
        df_0[continues_var].astype(float), 
        df_1[continues_var].astype(float),
        equal_var=False)
    output = [title, round(t,5), round(p,5), round((1-p)*100, 2), 
              round(df_0[continues_var].mean(),2),
              round(df_1[continues_var].mean(),2), 
              round(df_0[continues_var].mean(),2) - round(df_1[continues_var].mean(),2), 
              df_0.shape[0], df_1.shape[0]]
    print(significant(p, title), 
          '\nTest Statistic:', round(t,5), 
          '\np-value:', round(p,5),
          '\nConfidence:', round((1-p)*100, 2),
          '%\ncontinues_var Per User (Control): Â£', 
          round(df_0[continues_var].mean(),2),
          '\ncontinues_var Per User (Contender): Â£', 
          round(df_1[continues_var].mean(),2),
          '\ncontinues_var Difference:', 
          round((df_0[continues_var].mean() - df_1[continues_var].mean()),2),
          '\nSize (Control):', df_0.shape[0],
          '\nSize (Contender):', df_1.shape[0],   
          '\n--------------------------\n')
    return pd.DataFrame([output])

# read in data from DB
query = """SELECT variant_id, continues_var FROM table"""
conn = psycopg2.connect(
        host = host,
        password = password,
        port = port,
        user = user,
        database = database,
    )

cur = conn.cursor()
cur.execute(query)
df = pd.DataFrame(data = cur.fetchall(), 
                  columns = [column[0] for column in cur.description])
conn.close()
print('Number of row:,', df.shape[0])

# Calculate P-value and if it is significant
two_sided_t_test(df, 'title', 'continues_var', 'continues_var_name')
```
