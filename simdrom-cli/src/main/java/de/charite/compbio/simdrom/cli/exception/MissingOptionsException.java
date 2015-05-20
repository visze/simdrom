package de.charite.compbio.simdrom.cli.exception;

import java.util.List;

import org.apache.commons.lang3.StringUtils;

import de.charite.compbio.simdrom.cli.SIMdromSetting;

/**
 * Exception if a required option in {@link SIMdromSetting} is missing.
 * 
 * @author Max Schubach {@literal <max.schubach@charite.de>}
 *
 */
public class MissingOptionsException extends Exception {

	/**
	 * serial key for serialization
	 */
	private static final long serialVersionUID = -2694951897378018172L;

	/**
	 * Constructor Will build a message wit the two lists to notify the user.
	 * 
	 * @param usedOptions
	 *            List of necessary all options used.
	 * @param missingOptions
	 *            List of missing options.
	 */
	public MissingOptionsException(List<String> usedOptions, List<String> missingOptions) {
		super(buildMessage(usedOptions, missingOptions));

	}

	private static String buildMessage(List<String> usedOptions, List<String> missingOptions) {
		StringBuilder message = new StringBuilder("You cannot use the option(s) ");
		message.append(StringUtils.join(usedOptions, ", "));
		message.append(" without using option(s) ");
		message.append(StringUtils.join(usedOptions, ", "));
		return message.toString();
	}

}
