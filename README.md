# FRE6233 Final Exam

You have forked https://github.com/keithalewis/FRE6233-final 
to your GitHub repository and cloned that to your local machine.
Submit your solutions by pushing to your fork on GitHub and
issuing a pull request to have your work graded.

Written Exam

Write answers in this file or on paper.

Let (B_t) be standard Brownian motion.
The Bachelier model of stock price is F_t = f(1 + sigma B_t),
where f is the forward and sigma is the volatility

1. What is E[F_t] and Var(F_t)?

Solution:

E[F_t] = f since E[B_t] = 0, Var(F_t) = f^2 sigma^2 t since Var(B_t) = t.

Consider a European option that pays (F_T)^2 at time T.

2. Derive a formula for V_t = E_t[F_T^2] in terms of F_t.

Hint. Use B_t^2 - t is a martingale.

    E_t[(F_t)^2] = E_t[f^2(1 + 2 sigma B_T + sigma^2 B_T^2)]  
                 = f^2 E_t[1 + 2 sigma B_T + sigma^2 (B_T^2 - T + T)]  
                 = f^2 (1 + 2 sigma B_t + sigma^2 (B_t^2 - t + T)  
                 = f^2 (F_t)^2 + f^2 sigma^2 (T - t)  

Let Z be a standard normal random variable and N(z) = P(Z <= z) be its cumulative distribution function.

3. Show E[Z 1(Z <= a)] = -n(a) where n(z) = N'(z) is the standard normal density function.

Solution:

E[Z 1(Z <= a)] = int_-infty^a z e^{-z^2/2} dz/sqrt(2 pi)
               = -e^{z^2/2}/sqrt(2 pi) |_-infty^a
               = -n(a)

4. Prove the forward value of a European put option with strike k and expiration T is

    p = E[max{k - F_T, 0}] = (k - f) N(x) + f sigma sqrt(t) n(x)

where x = (k - f)/(f sigma sqrt(t)).

Hint. B_t/sqrt(t) is standard normal.

Solution:

    p = E[max{k - F_T, 0}]
      = E[(k - F_T)1(F_T <= k)]
      = k N(x) - E[F_t 1(F_t <= k)]
      = k N(x) - f N(x) + f sigma E[B_t 1(B_t <= (k - f)/f sigma)]
      = (k - f) N(x) + f sigma sqrt(t) n((k - f)/(f sigma sqrt(t)))

Since E[B_T 1(B_t/sqrt(T) <= x)] = sqrt(T) n(x).

5. Derive a closed form solution for the Bachelier delta dp/df.

Hint. Use dg(F_T)/df = g'(F_T) dF_T/df where g(x) = max{k - x, 0}.

Solution:

    dp/df = E[-1(F_T <= k)(1 + sigma B_T)]
          = -N(x) - sigma E[B_T 1(B_T/sqrt(T) <= x)]
          = -N(x) - sigma sqrt(T) n(x)

Coding Exam

Follow the instructions in xll_template.cpp.
