
package com.snuggy.nr.util;

import static java.lang.annotation.ElementType.*;

import java.lang.annotation.*;

@Retention(RetentionPolicy.RUNTIME)
@Target(value={CONSTRUCTOR, FIELD, LOCAL_VARIABLE, METHOD, PACKAGE, PARAMETER, TYPE})
public @interface Broken {
}
