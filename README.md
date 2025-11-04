# Be Bayesian my friend

Welcome to the short course **“Be Bayesian my friend”**, delivered at the **IV Congreso & XV Jornadas de Usuarios de R**.  

This course is designed for anyone interested in understanding and applying the principles of **Bayesian Inference** using **R** and modern computational tools.  
The goal is to bridge theory and practice, introducing participants to the Bayesian mindset and its implementation through some examples.

---

## 🧠 Course Overview

The workshop combines conceptual understanding with hands-on practice.  
It is divided into **two main parts**, to be covered in approximately **two hours**.

### **PART I — An Introduction to Bayesian Learning**
- History of Bayes’ theorem  
- Bayes’ theorem and its interpretation  
- Bayesian inference and posterior distribution  
- Credible intervals and predictive distributions  

### **PART II — Bayesian Computation and Mixed Models**
- From simple models to hierarchical structures  
- Bayesian computation: MCMC, Gibbs, and Metropolis–Hastings  
- Fitting models with `brms`  
- Practical applications.

---


## 💻 Software Requirements

To take full advantage of the course, please make sure you have the following software installed:

- [R (version 4.5.1 or later)](https://cran.r-project.org/) — the R 4.5.1 “Great Square Root” release or newer is recommended.  
- [RStudio](https://www.rstudio.com/products/rstudio/download/) — for a user-friendly coding interface.  
- [**JAGS**](https://sourceforge.net/projects/mcmc-jags/) — required for the `runjags` examples (used in the *S1-JAGS-heart_attack* practical exercise).  

> ⚠️ **Note:** `runjags` requires that the **JAGS** software is installed separately from the above link.

---

## 📦 R Packages

The following R packages are required for the theoretical and practical parts of the course:

```r
install.packages(pkgs = c(
  "ggplot2", "gridExtra", "dplyr", "readxl", "brms",
  "LearnBayes", "coda", "runjags"
))
```

Additional dependencies may be installed automatically when loading the main packages.

---

## 📁 Folder Structure

The repository is organized as follows:

```
Be_Bayesian_My_Friend/
│
├── Theory/
│   ├── Be Bayesian my Friend. Part I_ An Introduction to Bayesian Learning.html
│   ├── Be Bayesian my Friend. Part I_ An Introduction to Bayesian Learning.pdf
│   ├── Be Bayesian my Friend. Part II_ Bayesian Computation and Mixed Models.html
│   └── Be Bayesian my Friend. Part II_ Bayesian Computation and Mixed Models.pdf
│
├── Practical-lesson/
│   ├── S1-brms-heart_attack/
│   │   └── S1-brms-heart_attack.Rmd
│   │
│   └── S1-JAGS-heart_attack/
│       └── S1-JAGS-heart_attack.Rmd
│
└── README.md

```

- **Theory/** contains the conceptual slides and explanations (both in HTML and PDF).  
- **Practical-lesson/** includes applied exercises, such as a logistic regression example using the `brms` package.  

---

## 🎯 Learning Outcomes

By the end of this short course, participants will be able to:

- Understand the logic behind Bayesian inference and posterior reasoning.  
- Implement simple and hierarchical Bayesian models in R.  
- Use `brms` and other tools to fit Bayesian regression models.  
- Interpret posterior distributions, credible intervals, and predictions.  
- Recognize the value of Bayesian thinking in data analysis and decision making.

---

## 👨‍🏫 Instructor

**Joaquín Martínez-Minaya**  
Departamento de Estadística e Investigación Operativa Aplicadas y Calidad (DEIOAC)  
Universitat Politècnica de València (UPV)  

📧 **Email:** [jmarmin@eio.upv.es](mailto:jmarmin@eio.upv.es)  
🌐 **Website:** [https://github.com/jmartinez-minaya](https://github.com/jmartinez-minaya)


---

## 🔗 License

All materials are provided for educational purposes.  
Please acknowledge the author if you reuse or adapt any part of this course.
