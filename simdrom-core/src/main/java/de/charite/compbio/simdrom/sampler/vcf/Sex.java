/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

/**
 * @author max
 *
 */
public enum Sex {
	
	/**
	 * Sample a male genotype 
	 */
	MALE,
	/**
	 * Sample a female genotype
	 */
	FEMALE,
	/**
	 * No sex. gonosomes will be removed
	 */
	NONE,
	/**
	 * not known.
	 */
	UNKNOWN

}
