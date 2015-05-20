package de.charite.compbio.simdrom.cli.exception;

import java.util.List;

import org.apache.commons.cli.Option;
import org.apache.commons.lang3.StringUtils;

import de.charite.compbio.simdrom.cli.SIMdromSetting;

/**
 * Exception if a not allowed combination of {@link Option} is used in
 * {@link SIMdromSetting}.
 * 
 * @author Max Schubach {@literal <max.schubach@charite.de>}
 *
 */
public class NotAllowedCombinationOfOptionsException extends Exception {

	/**
	 * serial key for serialization
	 */
	private static final long serialVersionUID = -2694951897378018172L;

	/**
	 * Constructor builds a message with the given list.
	 * 
	 * @param values
	 *            List of options that are not allowed to combine
	 */
	public NotAllowedCombinationOfOptionsException(List<String> values) {
		super(buildMessage(values));

	}

	private static String buildMessage(List<String> values) {
		StringBuilder message = new StringBuilder("You cannot combine the options ");
		message.append(StringUtils.join(values, ", "));
		return message.toString();
	}

}
